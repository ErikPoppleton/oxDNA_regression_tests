#!/usr/bin/env python3
"""CPU/CUDA consistency of the capped backbone potential (max_backbone_force, max_backbone_rate).

The cap is computed by one class (Utilities/MaxBackboneForce) on the CPU, while each CUDA interaction evaluates the capped
potential in its own kernels from the parameters the class provides. This checks, for every interaction that supports
the cap (DNA, DNA2, DNA3, RNA2, NA, DNANM, RNANM, NANM and the ANNaMo plugin), on configurations where some backbone
bonds are beyond the cap:

  e0       energy of the GPU kernels (CUDA_print_energy) vs the host (CPU code) energy of the same configuration, with
           and without the cap: the difference with the cap must stay within the GPU-CPU difference without it
  ramp     max_backbone_rate on the GPU: with the step counter at T, the energy of max_backbone_force = F0 and
           max_backbone_rate = r must be the energy of a constant cap F0 + r * t, at the step t used by the host
  nve      constant-energy run with the cap: the drift and fluctuations of the conserved energy on the GPU must be
           comparable to the CPU ones (the GPU forces must be the derivatives of the capped potential)
  ramp_md  MD with a growing cap: GPU and host energies must agree at every printed step

  mbf_cuda_check.py prepare --cpu BIN [--plugins DIR]              (once: writes mbf/systems/*/e0.dat)
  mbf_cuda_check.py run --out DIR --cpu BIN --mixed BIN [--double BIN] [--cpu-plugins DIR] [--cuda-plugins DIR]
  mbf_cuda_check.py analyze --out DIR
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYS = os.path.join(REPO, 'mbf', 'systems')
sys.path.insert(0, REPO)
from regtest import oxio  # noqa: E402

COMMON = {'sim_type': 'MD', 'verlet_skin': 0.2, 'time_scale': 'linear', 'no_stdout_energy': 1, 'max_io': 100,
          'use_average_seq': 'true', 'energy_file': 'energy.dat', 'log_file': 'log.dat',
          'trajectory_file': 'trajectory.dat', 'lastconf_file': 'last_conf.dat', 'diff_coeff': 2.5,
          'newtonian_steps': 103}
CUDA = {'CUDA_list': 'verlet', 'CUDA_sort_every': 0, 'max_density_multiplier': 10}

# system: input keys, index of the nucleotide displaced to stretch its bonds (None: use the configuration as it is)
SYSTEMS = {
    'dsdna8': ({'T': '20C', 'salt_concentration': 0.5}, 3),
    'rna_duplex': ({'T': '20C', 'salt_concentration': 1.0}, 4),
    'na_mixed': ({'T': '25C', 'salt_concentration': 0.5}, 5),
    'ubq_dna': ({'T': '300K', 'salt_concentration': 1.0, 'parfile': 'anm.par'}, 80),
    'ubq_rna': ({'T': '300K', 'salt_concentration': 1.0, 'parfile': 'anm.par'}, 80),
    'ubq_hybrid': ({'T': '300K', 'salt_concentration': 1.0, 'parfile': 'anm.par'}, 80),
    'annamo': ({'T': 1.0, 'ANNAMO_annamo_version': 1, 'ANNAMO_tC': 30, 'ANNAMO_interaction_matrix_file': 'dHdS_matrix.dat',
                'DPS_3b_lambda': 1, 'generate_consider_bonded_interactions': 'true'}, None),
}
FILE_KEYS = ('parfile', 'ANNAMO_interaction_matrix_file')

# name, system, interaction, cap for e0/nve, CUDA precisions
CASES = [
    ('dsdna8_DNA', 'dsdna8', 'DNA', 1., ('mixed', 'double')),
    ('dsdna8_DNA2', 'dsdna8', 'DNA2', 1., ('mixed', 'double')),
    ('dsdna8_DNA3', 'dsdna8', 'DNA3', 1., ('mixed', 'double')),
    ('rna_duplex_RNA2', 'rna_duplex', 'RNA2', 1., ('mixed', 'double')),
    ('na_mixed_NA', 'na_mixed', 'NA', 1., ('mixed', 'double')),
    ('ubq_dna_DNANM', 'ubq_dna', 'DNANM', 1., ('mixed', 'double')),
    ('ubq_rna_RNANM', 'ubq_rna', 'RNANM', 1., ('mixed', 'double')),
    ('ubq_hybrid_NANM', 'ubq_hybrid', 'NANM', 1., ('mixed', 'double')),
    # the CUDA ANNaMo plugin does not compile in double precision
    ('annamo', 'annamo', 'ANNaMoInteraction', 5., ('mixed',)),
]
# ramp check: with the step counter at RAMP_T, the cap is F0 + rate * t = cap at t = RAMP_T and (1 + RAMP_FRAC) cap at RAMP_T + 1
RAMP_T, RAMP_FRAC = 2, 0.2   # F0 = 0.6 cap > 0


def ramp_params(cap):
    rate = RAMP_FRAC * cap
    return cap - rate * RAMP_T, rate


MD_F0, MD_RATE, MD_STEPS = 1., 0.002, 2000       # ramp_md: the cap grows from 1 to 5
NVE_STEPS, NVE_DT = 10000, 0.001
FLOOR = {'mixed': 2e-5, 'double': 1e-6}           # energies per particle


def write_input(path, keys):
    with open(path, 'w') as f:
        for k, v in keys.items():
            f.write(f'{k} = {v}\n')


def case_input(system, inter, backend, plugins):
    keys, _ = SYSTEMS[system]
    d = os.path.join(SYS, system)
    inp = {**COMMON, **keys, 'interaction_type': inter, 'topology': os.path.join(d, 'topology.top')}
    for k in FILE_KEYS:
        if k in inp:
            inp[k] = os.path.join(d, inp[k])
    if plugins:
        inp['plugin_search_path'] = plugins
    if backend == 'cpu':
        inp['backend'] = 'CPU'
    else:
        inp.update({'backend': 'CUDA', **CUDA, 'backend_precision': backend})
        # use_edge only in mixed precision; the ANNaMo CUDA plugin only implements the neighbour matrix
        inp['use_edge'] = 1 if backend == 'mixed' and system != 'annamo' else 0
        if inp['use_edge']:
            inp['edge_n_forces'] = 1
    return inp


def run(binary, rundir, inp, timeout=3600):
    os.makedirs(rundir, exist_ok=True)
    write_input(os.path.join(rundir, 'input'), inp)
    with open(os.path.join(rundir, 'stdout.txt'), 'w') as out:
        try:
            rc = subprocess.run([binary, 'input'], cwd=rundir, stdout=out, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    return rc


def conf_with_time(src, dst, t):
    lines = open(src).read().splitlines()
    lines[0] = f't = {t}'
    open(dst, 'w').write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------------------------------------------- prepare
def cmd_prepare(a):
    """Displace one nucleotide of each system so that some backbone bonds are beyond x_max for a cap of 2 but still
    within the FENE range, and check on the CPU that the cap changes the energy."""
    import shutil
    tmp = os.path.join(REPO, 'scratch', 'mbf_prepare')
    for system, (keys, idx) in SYSTEMS.items():
        d = os.path.join(SYS, system)
        conf = os.path.join(d, 'conf.dat')
        inter = next(c[2] for c in CASES if c[1] == system)
        cap = next(c[3] for c in CASES if c[1] == system)

        def energy(confpath, extra, inter=inter):
            rd = os.path.join(tmp, system)
            inp = {**case_input(system, inter, 'cpu', a.plugins), 'conf_file': confpath, 'steps': 0, 'dt': 0.001,
                   'thermostat': 'no', 'print_energy_every': 1, 'print_conf_interval': 1000, 'restart_step_counter': 1,
                   'refresh_vel': 1, 'seed': 1, **extra}
            if os.path.exists(os.path.join(rd, 'energy.dat')):
                os.remove(os.path.join(rd, 'energy.dat'))
            run(a.cpu, rd, inp)
            try:
                return oxio.read_energy(os.path.join(rd, 'energy.dat'))[0, 1]
            except Exception:
                return None

        if idx is None:
            shutil.copy(conf, os.path.join(d, 'e0.dat'))
            print(f'{system}: uncapped {energy(conf, {}):.6f}, capped {energy(conf, {"max_backbone_force": cap}):.6f}')
            continue
        inters = [c[2] for c in CASES if c[1] == system]
        lines = open(conf).read().splitlines()
        chosen = None
        for shift in (0.15, 0.12, 0.1, 0.08, 0.06, 0.04):
            for axis in (0, 1, 2):
                for sign in (1, -1):
                    v = lines[3 + idx].split()
                    v[axis] = str(float(v[axis]) + sign * shift)
                    new = lines[:3 + idx] + [' '.join(v)] + lines[4 + idx:]
                    path = os.path.join(d, 'e0.dat')
                    open(path, 'w').write('\n'.join(new) + '\n')
                    res = []
                    for inter in inters:
                        e_free, e_cap = energy(path, {}, inter), energy(path, {'max_backbone_force': cap}, inter)
                        res.append((inter, e_free, e_cap))
                    # every interaction must accept the configuration without the cap, and the cap must change the energy
                    if all(f is not None and c is not None and abs(f - c) > 1e-3 for _, f, c in res):
                        chosen = (f'{sign * shift:+g} along {"xyz"[axis]}', res)
                        break
                if chosen:
                    break
            if chosen:
                break
        if chosen is None:
            raise RuntimeError(f'{system}: no displacement gives a configuration within the FENE range with an active cap')
        print(f'{system}: nucleotide {idx} displaced by {chosen[0]}: ' + '; '.join(f'{i} uncapped {f:.6f} capped {c:.6f}' for i, f, c in chosen[1]))


# -------------------------------------------------------------------------------------------------------------------- run
def plan(a):
    bins = {'cpu': a.cpu, 'mixed': a.mixed, 'double': a.double}
    jobs = []
    for name, system, inter, cap, precisions in CASES:
        if a.only and name not in a.only:
            continue
        d = os.path.join(SYS, system)
        e0 = os.path.join(d, 'e0.dat')
        for backend in ('cpu',) + tuple(p for p in precisions if bins.get(p)):
            plugins = (a.cpu_plugins if backend == 'cpu' else a.cuda_plugins) if system == 'annamo' else None
            base = case_input(system, inter, backend, plugins)
            out = os.path.join(a.out, name, backend)
            one = {'steps': 0, 'dt': 0.001} if backend == 'cpu' else {'steps': 1, 'dt': 1e-6, 'CUDA_print_energy': 'true'}
            single = {**base, **one, 'thermostat': 'no', 'refresh_vel': 1, 'seed': 1, 'print_energy_every': 1,
                      'print_conf_interval': 1000}
            # e0
            jobs.append((bins[backend], os.path.join(out, 'e0_nocap'), {**single, 'conf_file': e0, 'restart_step_counter': 1}))
            jobs.append((bins[backend], os.path.join(out, 'e0_cap'), {**single, 'conf_file': e0, 'restart_step_counter': 1,
                                                                        'max_backbone_force': cap}))
            # ramp: the step counter starts at RAMP_T
            f0, rate = ramp_params(cap)
            for label, keys in (('ramp', {'max_backbone_force': f0, 'max_backbone_rate': rate}),
                                ('const_T', {'max_backbone_force': f0 + rate * RAMP_T}),
                                ('const_T1', {'max_backbone_force': f0 + rate * (RAMP_T + 1)})):
                rd = os.path.join(out, label)
                os.makedirs(rd, exist_ok=True)
                conf_with_time(e0, os.path.join(rd, 'start.dat'), RAMP_T)
                jobs.append((bins[backend], rd, {**single, 'conf_file': 'start.dat', 'restart_step_counter': 0, **keys}))
            # nve with a constant cap
            jobs.append((bins[backend], os.path.join(out, 'nve'), {**base, 'conf_file': e0, 'restart_step_counter': 1,
                                                                     'thermostat': 'no', 'refresh_vel': 1, 'seed': 7,
                                                                     'dt': NVE_DT, 'steps': NVE_STEPS,
                                                                     'print_energy_every': NVE_STEPS // 100,
                                                                     'print_conf_interval': NVE_STEPS,
                                                                     'max_backbone_force': cap}))
            # MD with a growing cap
            if backend != 'cpu':
                jobs.append((bins[backend], os.path.join(out, 'ramp_md'), {**base, 'conf_file': e0, 'restart_step_counter': 1,
                                                                             'thermostat': 'john', 'refresh_vel': 1, 'seed': 9,
                                                                             'dt': 0.001, 'steps': MD_STEPS,
                                                                             'print_energy_every': MD_STEPS // 20,
                                                                             'print_conf_interval': MD_STEPS,
                                                                             'CUDA_print_energy': 'true',
                                                                             'max_backbone_force': MD_F0,
                                                                             'max_backbone_rate': MD_RATE}))
    return jobs


def cmd_run(a):
    a.out = os.path.abspath(a.out)
    jobs = plan(a)
    print(f'{len(jobs)} runs', flush=True)
    failed = 0
    for i, (binary, rd, inp) in enumerate(jobs, 1):
        rc = run(binary, rd, inp)
        if rc != 0:
            failed += 1
            err = ''
            log = os.path.join(rd, 'log.dat')
            if os.path.exists(log):
                errs = [l.strip() for l in open(log) if 'ERROR' in l]
                err = errs[-1] if errs else ''
            print(f'FAILED ({rc}) {rd}: {err}', flush=True)
        if i % 20 == 0:
            print(f'{i}/{len(jobs)} done', flush=True)
    print(f'{failed} failed runs')
    return 1 if failed else 0


# ---------------------------------------------------------------------------------------------------------------- analyze
def energies(rd, gpu=False):
    path = os.path.join(rd, 'energy.dat')
    if not os.path.exists(path):
        return None
    en, g = oxio.read_energy(path, with_gpu=True)
    if len(en) == 0:
        return None
    return (en, np.array(g)) if gpu else en


def nve_stats(rd):
    en = energies(rd)
    if en is None or len(en) < 3:
        return None
    t, E = en[:, 0], en[:, 1] + en[:, 2]
    slope, icpt = np.polyfit(t, E, 1)
    return {'drift': float(slope * (t[-1] - t[0])), 'fluct': float((E - slope * t - icpt).std()),
            'finite': bool(np.all(np.isfinite(E)))}


def cmd_analyze(a):
    a.out = os.path.abspath(a.out)
    rows, results = [], {}

    def add(case, backend, check, ok, msg):
        status = 'PASS' if ok else ('FAIL' if ok is not None else 'ERROR')
        rows.append((case, backend, check, status, msg))
        results.setdefault(case, {}).setdefault(backend, {})[check] = {'status': status, 'message': msg}

    for name, system, inter, cap, precisions in CASES:
        cdir = os.path.join(a.out, name)
        if not os.path.isdir(cdir):
            continue
        cpu = os.path.join(cdir, 'cpu')
        cpu_nve = nve_stats(os.path.join(cpu, 'nve'))
        e = {k: energies(os.path.join(cpu, k)) for k in ('e0_nocap', 'e0_cap', 'ramp', 'const_T')}
        if all(v is not None for v in e.values()):
            ok = e['ramp'][0, 1] == e['const_T'][0, 1] and e['e0_cap'][0, 1] != e['e0_nocap'][0, 1]
            add(name, 'cpu', 'ramp', ok, f'E(cap active) {e["e0_cap"][0, 1]:.6f} vs uncapped {e["e0_nocap"][0, 1]:.6f}; '
                f'ramp at t = {RAMP_T}: {e["ramp"][0, 1]:.6f}, constant {cap:g}: {e["const_T"][0, 1]:.6f}')
        else:
            add(name, 'cpu', 'ramp', None, 'missing CPU energies')
        for prec in precisions:
            pdir = os.path.join(cdir, prec)
            if not os.path.isdir(pdir):
                continue
            floor = FLOOR[prec]
            # e0: GPU vs host on the same configuration (row 1 of a one-step run)
            g = {k: energies(os.path.join(pdir, k), gpu=True) for k in ('e0_nocap', 'e0_cap', 'ramp', 'const_T', 'const_T1')}
            if any(v is None or len(v[0]) < 2 or not np.isfinite(v[1][1]) for v in g.values()):
                add(name, prec, 'e0', None, 'missing GPU energies (CUDA_print_energy)')
                continue
            d_free = g['e0_nocap'][1][1] - g['e0_nocap'][0][1, 1]
            d_cap = g['e0_cap'][1][1] - g['e0_cap'][0][1, 1]
            tol = 3 * abs(d_free) + floor
            add(name, prec, 'e0', abs(d_cap) <= tol,
                f'GPU - host: capped {d_cap:+.2e}, uncapped {d_free:+.2e} (tolerance {tol:.1e}); '
                f'capped host {g["e0_cap"][0][1, 1]:.6f} (CPU backend {e["e0_cap"][0, 1]:.6f})' if e['e0_cap'] is not None else '')
            # ramp: host and GPU must use the cap of the same step
            host = {k: g[k][0][1, 1] for k in ('ramp', 'const_T', 'const_T1')}
            gpu = {k: g[k][1][1] for k in ('ramp', 'const_T', 'const_T1')}
            k_host = min(('const_T', 'const_T1'), key=lambda k: abs(host[k] - host['ramp']))
            sep = abs(gpu['const_T'] - gpu['const_T1'])
            ok = abs(gpu['ramp'] - gpu[k_host]) <= max(floor, 1e-3 * sep) and sep > 10 * floor and abs(host['ramp'] - host[k_host]) < 1e-9
            add(name, prec, 'ramp', ok,
                f'host uses the cap of step {RAMP_T if k_host == "const_T" else RAMP_T + 1}; GPU ramp {gpu["ramp"]:.6f}, '
                f'GPU constant caps {gpu["const_T"]:.6f} (t = {RAMP_T}) / {gpu["const_T1"]:.6f} (t = {RAMP_T + 1})')
            # nve
            gn = nve_stats(os.path.join(pdir, 'nve'))
            if gn is None or cpu_nve is None:
                add(name, prec, 'nve', None, 'missing NVE energies')
            else:
                nfloor = 3e-4 if prec == 'mixed' else 1e-4
                ok = gn['finite'] and abs(gn['drift']) <= 3 * abs(cpu_nve['drift']) + nfloor and gn['fluct'] <= 3 * cpu_nve['fluct'] + nfloor
                add(name, prec, 'nve', ok, f'drift GPU {gn["drift"]:+.2e} CPU {cpu_nve["drift"]:+.2e}, fluct GPU {gn["fluct"]:.1e} '
                    f'CPU {cpu_nve["fluct"]:.1e} (per particle, {NVE_STEPS * NVE_DT:g} time units)')
            # ramp_md
            r = energies(os.path.join(pdir, 'ramp_md'), gpu=True)
            if r is None:
                add(name, prec, 'ramp_md', None, 'missing energies')
            else:
                en, gg = r
                ok_rows = np.isfinite(gg)
                diff = np.abs(gg[ok_rows] - en[ok_rows, 1])
                ok = len(diff) > 0 and diff.max() <= tol
                add(name, prec, 'ramp_md', ok, f'{len(diff)} steps, max |GPU - host| {diff.max() if len(diff) else float("nan"):.2e} '
                    f'(tolerance {tol:.1e}) while the cap grows from {MD_F0:g} to {MD_F0 + MD_RATE * MD_STEPS:g}')

    lines = ['# max_backbone_force CPU/CUDA consistency', '',
             '| case | backend | check | status | details |', '|---|---|---|---|---|']
    lines += [f'| {c} | {b} | {k} | {s} | {m} |' for c, b, k, s, m in rows]
    n_bad = sum(r[3] != 'PASS' for r in rows)
    lines += ['', f'{len(rows) - n_bad} passed, {n_bad} failed or missing']
    text = '\n'.join(lines) + '\n'
    open(os.path.join(a.out, 'mbf_cuda_check.md'), 'w').write(text)
    json.dump(results, open(os.path.join(a.out, 'mbf_cuda_check.json'), 'w'), indent=1)
    print(text)
    return 1 if n_bad else 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    pr = sub.add_parser('prepare')
    pr.add_argument('--cpu', required=True)
    pr.add_argument('--plugins', help='folder with ANNaMoInteraction.so')
    r = sub.add_parser('run')
    r.add_argument('--out', required=True)
    r.add_argument('--cpu', required=True)
    r.add_argument('--mixed')
    r.add_argument('--double')
    r.add_argument('--cpu-plugins', help='folder with the CPU ANNaMoInteraction.so')
    r.add_argument('--cuda-plugins', help='folder with ANNaMoInteraction.so and CUDAANNaMoInteraction.so of the CUDA build')
    r.add_argument('--only', nargs='*')
    an = sub.add_parser('analyze')
    an.add_argument('--out', required=True)
    a = p.parse_args()
    for k in ('cpu', 'mixed', 'double'):
        if getattr(a, k, None):
            setattr(a, k, os.path.abspath(getattr(a, k)))
    return {'prepare': cmd_prepare, 'run': cmd_run, 'analyze': cmd_analyze}[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
