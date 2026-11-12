#!/usr/bin/env python3
"""Equipartition check of the oxDNA MD backends (CPU vs CUDA, current branch vs master).

The CUDA regression suite found the GPU runs of protein systems a few percent colder than the CPU
ones (both the kinetic and the potential energy low by the same fraction, NA-only systems fine).
With the Brownian ("john") thermostat the kinetic energy must sample kT/2 per degree of freedom
whatever the forces are, as long as the dynamics between refreshes is Hamiltonian, so this script
measures the temperature directly from the velocities and angular momenta written to the
trajectory, split by particle type (amino acids / nucleotides) and by translation / rotation:

    T_trans = <v^2>/3      T_rot = <L^2>/3      (mass = inertia moment = 1)

Runs (each with many independent replicas, see RUNS below):
  cln025 DNANM / DNACT     10-residue protein, the coldest system in the suite
  cln025 DNACT, dt/2       a time-step artefact shrinks ~4x, a plain bug does not
  cln025 DNACT, pt = 1     every velocity refreshed every step: tests thermostat + output only
  ubq_dna DNANM            protein + DNA: are the nucleotides cold too?
  na_dna_duplex NA         NA-only control
on the CPU, CUDA mixed and CUDA double builds of the branch, and CUDA DNANM of master (if given).

  equipartition.py run --out DIR --branch-cpu BIN --branch-mixed BIN [--branch-double BIN] [--master-mixed BIN]
  equipartition.py analyze --out DIR
"""
import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYS = os.path.join(REPO, 'systems')

BASE = {
    'sim_type': 'MD', 'verlet_skin': 0.2, 'thermostat': 'john', 'refresh_vel': 1, 'time_scale': 'linear',
    'restart_step_counter': 1, 'no_stdout_energy': 1, 'max_io': 100,
    'trajectory_file': 'trajectory.dat', 'lastconf_file': 'last_conf.dat', 'energy_file': 'energy.dat',
    'log_file': 'log.dat',
}
CUDA = {'CUDA_list': 'verlet', 'CUDA_sort_every': 0, 'max_density_multiplier': 10}

PROTEIN = {'T': '300K', 'dt': 0.002, 'newtonian_steps': 51, 'diff_coeff': 1, 'salt_concentration': 1.0}

# system: input keys, protein / nucleic particle ranges
SYSTEMS = {
    'cln025': ({'conf_file': 'equilibrated.dat', 'topology': 'mixed.top', **PROTEIN}, [0, 10], None),
    'ubq_dna': ({'conf_file': 'equilibrated.dat', 'topology': 'hybrid.top', 'external_forces': 1,
                 'external_forces_file': 'forces.txt', **PROTEIN}, [0, 76], [76, 100]),
    'na_dna_duplex': ({'conf_file': 'conf.dat', 'topology': 'topology.top', 'interaction_type': 'NA',
                       'salt_concentration': 0.5, 'use_average_seq': 'true', 'T': '25C', 'dt': 0.002,
                       'diff_coeff': 2.5, 'newtonian_steps': 103}, None, [0, 22]),
}
FILE_KEYS = ('conf_file', 'topology', 'parfile', 'external_forces_file')

DNANM = {'interaction_type': 'DNANM', 'parfile': 'anm.par'}
DNACT = {'interaction_type': 'DNACT', 'parfile': 'anmt.par'}

# name, system, extra input keys, steps, conf/energy interval, replicas, builds
ALL = ('branch_cpu', 'branch_mixed', 'branch_double')
RUNS = [
    ('cln025_DNANM', 'cln025', DNANM, 300000, 500, 20, ALL + ('master_mixed',)),
    ('cln025_DNACT', 'cln025', DNACT, 300000, 500, 20, ALL),
    ('cln025_DNACT_dt0.001', 'cln025', {**DNACT, 'dt': 0.001, 'newtonian_steps': 102}, 600000, 1000, 20, ALL),
    ('cln025_DNACT_pt1', 'cln025', {**DNACT, 'newtonian_steps': 1, 'pt': 1.0}, 50000, 100, 5, ALL),
    ('ubq_dna_DNANM', 'ubq_dna', DNANM, 80000, 400, 12, ('branch_cpu', 'branch_mixed', 'master_mixed')),
    ('na_dna_duplex_NA', 'na_dna_duplex', {}, 100000, 500, 12, ('branch_cpu', 'branch_mixed')),
]
EQUIL = 0.2   # fraction of each run discarded


def temperature(t):
    t = str(t).strip()
    if t[-1] in 'kK':
        return float(t[:-1]) / 3000.
    if t[-1] in 'cC':
        return (float(t[:-1]) + 273.15) / 3000.
    return float(t)


def make_input(rundir, system, extra, steps, every, seed, build):
    keys, _, _ = SYSTEMS[system]
    inp = {**BASE, **keys, **extra, 'steps': steps, 'print_energy_every': every, 'print_conf_interval': every,
           'seed': seed}
    for k in FILE_KEYS:
        if k in inp:
            inp[k] = os.path.join(SYS, system, inp[k])
    if build == 'branch_cpu':
        inp.update({'backend': 'CPU', 'backend_precision': 'double'})
    else:
        inp.update({'backend': 'CUDA', **CUDA})
        if build.endswith('double'):
            inp.update({'backend_precision': 'double', 'use_edge': 0})
        else:
            # the branch suite ran mixed with use_edge = 1; master's DNANM is run with the default kernels
            inp.update({'backend_precision': 'mixed', 'use_edge': 0 if build.startswith('master') else 1})
            if not build.startswith('master'):
                inp['edge_n_forces'] = 1
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(rundir, 'input'), 'w') as f:
        for k, v in inp.items():
            f.write(f'{k} = {v}\n')


def run_one(binary, rundir, timeout):
    with open(os.path.join(rundir, 'stdout.txt'), 'w') as out:
        try:
            rc = subprocess.run([binary, 'input'], cwd=rundir, stdout=out, stderr=subprocess.STDOUT,
                                timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
        except OSError as ex:
            rc = str(ex)
    return rundir, rc


def cmd_run(a):
    bins = {'branch_cpu': a.branch_cpu, 'branch_mixed': a.branch_mixed, 'branch_double': a.branch_double,
            'master_mixed': a.master_mixed}
    cpu_jobs, gpu_jobs = [], []
    for name, system, extra, steps, every, reps, builds in RUNS:
        if a.only and name not in a.only:
            continue
        for b in builds:
            if not bins.get(b):
                continue
            for r in range(reps):
                rundir = os.path.join(a.out, name, b, f'rep{r}')
                if os.path.isfile(os.path.join(rundir, 'last_conf.dat')) and not a.force:
                    continue
                make_input(rundir, system, extra, max(int(steps * a.scale), every * 10), every, 1000 + 37 * r + (5000 if b != 'branch_cpu' else 0), b)
                (cpu_jobs if b == 'branch_cpu' else gpu_jobs).append((bins[b], rundir))
    print(f'{len(cpu_jobs)} CPU and {len(gpu_jobs)} GPU runs', flush=True)
    failed = []
    with cf.ThreadPoolExecutor(a.jobs) as cpu, cf.ThreadPoolExecutor(a.gpu_jobs) as gpu:
        futs = [cpu.submit(run_one, b, d, a.timeout) for b, d in cpu_jobs]
        futs += [gpu.submit(run_one, b, d, a.timeout) for b, d in gpu_jobs]
        for i, f in enumerate(cf.as_completed(futs), 1):
            d, rc = f.result()
            if rc != 0:
                failed.append(d)
                print(f'FAILED ({rc}): {d}', flush=True)
            if i % 50 == 0:
                print(f'{i}/{len(futs)} done', flush=True)
    print(f'{len(failed)} failed runs')
    return 1 if failed else 0


def read_traj(path):
    """(n_frames, N, 6): velocities and angular momenta of every frame."""
    frames, cur = [], []
    with open(path) as f:
        for line in f:
            c = line[:3]
            if c == 't =':
                if cur:
                    frames.append(cur)
                cur = []
            elif c in ('b =', 'E ='):
                continue
            elif line.strip():
                cur.append([float(x) for x in line.split()[9:15]])
    if cur:
        frames.append(cur)
    n = min(len(fr) for fr in frames)
    return np.array([fr for fr in frames if len(fr) == n])


def analyze_rep(rundir, T):
    x = read_traj(os.path.join(rundir, 'trajectory.dat'))
    x = x[int(len(x) * EQUIL):]
    v2 = (x[:, :, 0:3] ** 2).mean(axis=(0, 2)) / T      # per particle <v_a^2>/T: 1 at equipartition
    l2 = (x[:, :, 3:6] ** 2).mean(axis=(0, 2)) / T
    e = np.loadtxt(os.path.join(rundir, 'energy.dat'), ndmin=2)
    e = e[int(len(e) * EQUIL):]
    return v2, l2, float(e[:, 1].mean()), float(e[:, 2].mean())


def mean_se(vals):
    vals = np.asarray(vals, float)
    if len(vals) < 2:
        return float(vals.mean()) if len(vals) else float('nan'), float('nan')
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(len(vals)))


def cmd_analyze(a):
    results, lines = {}, []
    lines += ['# Equipartition check', '',
              'T_trans = <v^2>/3 and T_rot = <L^2>/3 over the frames after the first '
              f'{int(EQUIL * 100)}% of each replica, divided by the thermostat T (1 = equipartition); '
              '+- is the standard error over replicas, z the deviation from 1 in units of it. '
              'Rotational temperatures of 0 mean the particles do not rotate (e.g. DNANM amino acids). '
              'U and K are the per-particle potential and kinetic energies of energy.dat.', '']
    lines += ['| run | build | reps | group | T_trans/T | z | T_rot/T | z | U | K |', '|---|---|---|---|---|---|---|---|---|---|']
    for name, system, extra, steps, every, reps, builds in RUNS:
        keys, prot, nuc = SYSTEMS[system]
        T = temperature({**keys, **extra}['T'])
        for b in builds:
            bdir = os.path.join(a.out, name, b)
            if not os.path.isdir(bdir):
                continue
            per = []
            for r in sorted(os.listdir(bdir)):
                d = os.path.join(bdir, r)
                try:
                    per.append(analyze_rep(d, T))
                except Exception as ex:  # noqa: BLE001  (a failed run should not stop the analysis)
                    print(f'skipping {d}: {ex}', file=sys.stderr)
            if not per:
                lines.append(f'| {name} | {b} | 0 | - | no data | | | | | |')
                continue
            v2 = np.array([p[0] for p in per])
            l2 = np.array([p[1] for p in per])
            U = mean_se([p[2] for p in per])
            K = mean_se([p[3] for p in per])
            np.savez(os.path.join(bdir, 'per_particle.npz'), v2=v2, l2=l2, U=[p[2] for p in per],
                     K=[p[3] for p in per])
            groups = [('protein', prot), ('nucleic', nuc)]
            groups = [(g, rng) for g, rng in groups if rng]
            res = {'replicas': len(per), 'U': U, 'K': K}
            for g, (i0, i1) in groups:
                tt = mean_se(v2[:, i0:i1].mean(axis=1))
                tr = mean_se(l2[:, i0:i1].mean(axis=1))
                res[g] = {'T_trans': tt, 'T_rot': tr}
                zt = (tt[0] - 1) / tt[1] if tt[1] > 0 else float('nan')
                zr = (tr[0] - 1) / tr[1] if tr[1] > 0 and tr[0] > 1e-6 else float('nan')
                lines.append(f'| {name} | {b} | {len(per)} | {g} | {tt[0]:.4f}+-{tt[1]:.4f} | {zt:+.1f} | '
                             f'{tr[0]:.4f}+-{tr[1]:.4f} | {zr:+.1f} | {U[0]:+.5f}+-{U[1]:.5f} | {K[0]:.5f}+-{K[1]:.5f} |')
            results.setdefault(name, {})[b] = res
    with open(os.path.join(a.out, 'equipartition.json'), 'w') as f:
        json.dump(results, f, indent=1)
    text = '\n'.join(lines) + '\n'
    with open(os.path.join(a.out, 'equipartition.md'), 'w') as f:
        f.write(text)
    print(text)
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--out', required=True)
    r.add_argument('--branch-cpu', required=True)
    r.add_argument('--branch-mixed', required=True)
    r.add_argument('--branch-double')
    r.add_argument('--master-mixed')
    r.add_argument('--jobs', type=int, default=8)
    r.add_argument('--gpu-jobs', type=int, default=4)
    r.add_argument('--timeout', type=float, default=3600)
    r.add_argument('--only', nargs='*', help='run names (see RUNS)')
    r.add_argument('--scale', type=float, default=1.0, help='scale the number of steps (for testing)')
    r.add_argument('--force', action='store_true', help='rerun replicas that already finished')
    an = sub.add_parser('analyze')
    an.add_argument('--out', required=True)
    a = p.parse_args()
    a.out = os.path.abspath(a.out)
    return cmd_run(a) if a.cmd == 'run' else cmd_analyze(a)


if __name__ == '__main__':
    sys.exit(main())
