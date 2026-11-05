#!/usr/bin/env python3
"""Regression tests for the oxDNA protein (ANM) models and for oxNA.

Each case (systems/<system>/config.json -> "cases") is run on a reference side ("old") and a test side
("new") and compared through up to four stages:

  energy0  potential energy of a perturbed configuration (steps = 0); also checked against an
           independent numpy evaluation of the protein terms (tools/anm_reference.py). On CUDA the
           energy computed on the GPU (CUDA_print_energy) is also compared with the host one.
  det      short MD run with the same seed and thermostat: trajectories must be (near) identical
           (only when both sides use the same backend)
  stat     independent replicas: mean potential energy and per-residue RMSF within statistical error
  nve      constant-energy run (thermostat = no): drift and fluctuations of the conserved energy

Two modes:
  original vs current (default): --old is the original anm-oxdna binary, --new the current oxDNA.
      Cases with "old_code": "current" run their "old" inputs with the --new binary instead (e.g. NANM
      vs DNANM), cases with "old_code": "stored" compare with values stored in systems/<s>/stored.json.
  same code (--same-code): both sides use the "new" inputs of every case; typically --old is a CPU
      build and --new a CUDA build of the same oxDNA (--backend CUDA).

Examples:
  ./run_tests.py --old OLD/bin/oxDNA --new NEW/bin/oxDNA --quick
  ./run_tests.py --same-code --old CPU/bin/oxDNA --new CUDA/bin/oxDNA --backend CUDA --cuda-precision mixed
"""
import argparse
import concurrent.futures as cf
import fnmatch
import json
import os
import shutil
import subprocess
import sys
import threading
import time

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tools'))
from regtest import oxio, analysis  # noqa: E402
import anm_reference  # noqa: E402

BASE_INPUT = {
    'sim_type': 'MD',
    'T': '300K',
    'dt': 0.002,
    'verlet_skin': 0.2,
    'thermostat': 'john',
    'newtonian_steps': 51,
    'diff_coeff': 1,
    'refresh_vel': 1,
    'salt_concentration': 1.0,
    'time_scale': 'linear',
    'restart_step_counter': 1,
    'no_stdout_energy': 1,
    'max_io': 100,
}
CUDA_COMMON = {'backend': 'CUDA', 'CUDA_list': 'verlet', 'CUDA_sort_every': 0, 'max_density_multiplier': 10}
BACKEND_INPUT = {
    'CPU': {'backend': 'CPU', 'backend_precision': 'double'},
    # the ANM interactions require use_edge on CUDA; use_edge is not available in double precision
    'CUDA': dict(CUDA_COMMON, backend_precision='mixed', use_edge=1, edge_n_forces=1),
    'CUDA_float': dict(CUDA_COMMON, backend_precision='float', use_edge=1, edge_n_forces=1),
    'CUDA_double': dict(CUDA_COMMON, backend_precision='double', use_edge=0),
}
FILE_KEYS = ('topology', 'parfile', 'conf_file', 'external_forces_file', 'massfile')
ALL_STAGES = ['energy0', 'det', 'stat', 'nve']

PASS, FAIL, XFAIL, XPASS, SKIP, ERROR = 'PASS', 'FAIL', 'XFAIL', 'XPASS', 'SKIP', 'ERROR'


# ----------------------------------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------------------------------

def load_cases(args):
    cases = []
    sysdir = os.path.join(ROOT, 'systems')
    for name in sorted(os.listdir(sysdir)):
        cfg_path = os.path.join(sysdir, name, 'config.json')
        if not os.path.isfile(cfg_path):
            continue
        with open(cfg_path) as f:
            cfg = json.load(f)
        for cname, case in cfg['cases'].items():
            full = f'{name}/{cname}'
            if args.only and not any(fnmatch.fnmatch(full, pat) or fnmatch.fnmatch(name, pat) for pat in args.only):
                continue
            if args.quick and not case.get('quick', False):
                continue
            if args.smoke and not case.get('smoke', False):
                continue
            if args.mode not in case.get('modes', ['original', 'same_code']):
                continue
            if args.backend == 'CPU' and case.get('cuda_only', False):
                continue
            c = dict(case)
            c['system'] = name
            c['name'] = full
            c['dir'] = os.path.join(sysdir, name)
            c['sys'] = cfg
            cases.append(c)
    return cases


def stage_settings(case, stage, args):
    s = dict(case['sys'].get('stages', {}).get(stage, {}))
    s.update(case.get('stages', {}).get(stage, {}))
    for sub in ('quick', 'smoke'):
        o = s.pop(sub, {})
        if getattr(args, sub):
            s.update(o)
    return s


def enabled_stages(case, args):
    st = case.get(f'stages_enabled@{args.mode}', case.get('stages_enabled', ALL_STAGES))
    return [x for x in args.stages if x in st]


def expectation(case, stage, args):
    """Return the reason string if this check is expected to fail with this backend, else None."""
    tags = [args.btag] + (['CUDA'] if args.btag.startswith('CUDA') and args.btag != 'CUDA' else [])
    keys = []
    for t in tags:
        keys += [f'{stage}@{t}', f'*@{t}']
    keys += [f'{stage}@{args.mode}', f'*@{args.mode}', stage, '*']
    for key in keys:
        if key in case.get('xfail', {}):
            return case['xfail'][key]
    return None


def _xfail_applies(x, chk):
    """An expected failure is either a reason string or {"reason": ..., "match": ...}: with "match",
    only a failure whose message contains that text counts as expected."""
    if isinstance(x, dict):
        return x.get('match') is None or x['match'] in chk['message'], x['reason']
    return True, x


def side_plan(case, side, args):
    """(binary, backend tag, list of input-set names) for one side, or None if the side is not run."""
    test_sets = ['new', f'new@{args.btag}'] + (['new@CUDA'] if args.btag.startswith('CUDA_') else [])
    if side == 'new':
        return args.new, args.btag, test_sets
    if args.mode == 'same_code':
        # identical inputs on both sides (including the CUDA-specific substitutions), only the backend differs
        return args.old, args.old_btag, test_sets
    code = case.get('old_code', 'original')
    if code == 'stored':
        return None
    if code == 'current':
        return args.new, 'CPU', ['old', 'old@CPU']
    return args.old, args.old_btag, ['old', f'old@{args.old_btag}']


def make_input(case, input_sets, btag, stage_overrides, rundir):
    inp = dict(BASE_INPUT)
    inp.update(BACKEND_INPUT[btag])
    inp.update(case['sys'].get('common', {}))
    inp.update(case.get('common', {}))
    for name in input_sets:
        inp.update(case.get(name, {}))
    inp.update(stage_overrides)
    for k in FILE_KEYS:
        if k in inp and not os.path.isabs(str(inp[k])):
            inp[k] = os.path.join(case['dir'], inp[k])
    inp.update({'trajectory_file': 'trajectory.dat', 'lastconf_file': 'last_conf.dat',
                'energy_file': 'energy.dat', 'log_file': 'log.dat'})
    os.makedirs(rundir, exist_ok=True)
    blocks = inp.pop('_data_output', [])
    with open(os.path.join(rundir, 'input'), 'w') as f:
        for k, v in inp.items():
            f.write(f'{k} = {v}\n')
        for i, b in enumerate(blocks, 1):
            f.write(f'data_output_{i} = {{\n{b}\n}}\n')
    return inp


def run_oxdna(binary, rundir, timeout, gpu_sem=None):
    if gpu_sem is not None:
        gpu_sem.acquire()
    t0 = time.time()
    try:
        with open(os.path.join(rundir, 'stdout.txt'), 'w') as out:
            try:
                p = subprocess.run([binary, 'input'], cwd=rundir, stdout=out, stderr=subprocess.STDOUT,
                                   timeout=timeout)
                rc = p.returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
    finally:
        if gpu_sem is not None:
            gpu_sem.release()
    elapsed = time.time() - t0
    txt = ''
    for name in ('stdout.txt', 'log.dat'):  # errors go to log_file when it is set
        if os.path.isfile(os.path.join(rundir, name)):
            with open(os.path.join(rundir, name), errors='replace') as f:
                txt += f.read()
    err = None
    if rc != 0 or 'ERROR' in txt:
        lines = [l for l in txt.splitlines() if 'ERROR' in l]
        err = lines[-1].strip() if lines else f'exit code {rc}'
    return {'rc': rc, 'time': elapsed, 'error': err}


# ----------------------------------------------------------------------------------------------------
# stored references (cases with "old_code": "stored")
# ----------------------------------------------------------------------------------------------------

def _stored_path(case):
    return os.path.join(case['dir'], 'stored.json')


def stored_value(case, key):
    p = _stored_path(case)
    if not os.path.isfile(p):
        return None
    with open(p) as f:
        return json.load(f).get(case['name'].split('/', 1)[1], {}).get(key)


def update_stored(case, key, value, info):
    p = _stored_path(case)
    data = {}
    if os.path.isfile(p):
        with open(p) as f:
            data = json.load(f)
    entry = data.setdefault(case['name'].split('/', 1)[1], {})
    entry[key] = value
    entry['_generated_with'] = info
    with open(p, 'w') as f:
        json.dump(data, f, indent=1, default=float)


# ----------------------------------------------------------------------------------------------------
# stages
# ----------------------------------------------------------------------------------------------------

def _has_ext_forces(case, input_sets):
    inp = dict(case['sys'].get('common', {}))
    inp.update(case.get('common', {}))
    for n in input_sets:
        inp.update(case.get(n, {}))
    return str(inp.get('external_forces', 0)).lower() in ('1', 'true', 'yes')


def stage_runs(case, stage, args):
    """Return the list of (side, rundir, input-overrides) for a stage."""
    s = stage_settings(case, stage, args)
    base = os.path.join(args.out, case['name'], stage)
    sides = [x for x in ('old', 'new') if side_plan(case, x, args) is not None]
    runs = []
    if stage == 'energy0':
        # external forces off: only the interaction energy is compared
        ov = {'steps': 0, 'print_energy_every': 1, 'print_conf_interval': 1000, 'seed': 1, 'external_forces': 0,
              'conf_file': case['sys'].get('energy0_conf', case['sys']['common']['conf_file'])}
        for side in sides:
            o0 = dict(ov)
            if side_plan(case, side, args)[1].startswith('CUDA'):
                # one step: the line at t = dt carries both the host energy and the GPU one
                o0.update({'steps': 1, 'CUDA_print_energy': 'true', 'thermostat': 'no', 'dt': 1e-6})
            runs.append((side, os.path.join(base, side), o0))
            for vname, vov in case.get('energy0_variants', {}).items():
                o = dict(o0)
                o.update(vov)
                runs.append((side, os.path.join(base, f'{side}_{vname}'), o))
    elif stage == 'det':
        n = s['steps']
        ov = {'steps': n, 'print_energy_every': s.get('print_every', n // 5),
              'print_conf_interval': s.get('print_every', n // 5), 'seed': 1}
        for side in sides:
            runs.append((side, os.path.join(base, side), ov))
    elif stage == 'stat':
        # stored references are generated with the seeds of the reference side, so that later runs of the
        # test side (seeds 2000+) are statistically independent of them
        new_seed = 1000 if (args.update_stored and case.get('old_code') == 'stored') else 2000
        for side, seed in (('old', 1000), ('new', new_seed)):
            if side not in sides:
                continue
            for vname, vov in [('', {})] + list(case.get('stat_variants', {}).items()):
                for rep in range(s.get('replicas', 1)):
                    ov = {'steps': s['steps'], 'print_energy_every': s['print_energy_every'],
                          'print_conf_interval': s['print_conf_interval'], 'seed': seed + rep}
                    if 'stat_conf' in case['sys']:
                        ov['conf_file'] = case['sys']['stat_conf']
                    ov.update(vov)
                    runs.append((side, os.path.join(base, side + (f'_{vname}' if vname else ''), f'rep{rep}'), ov))
    elif stage == 'nve':
        n = s.get('steps', 20000)
        every = s.get('print_every', max(1, n // 200))
        for side in sides:
            ov = {'steps': n, 'thermostat': 'no', 'refresh_vel': 1, 'seed': 7, 'dt': s.get('dt', 0.001),
                  'print_energy_every': every, 'print_conf_interval': n}
            if 'stat_conf' in case['sys']:
                ov['conf_file'] = case['sys']['stat_conf']
            if _has_ext_forces(case, side_plan(case, side, args)[2]):
                ov['_data_output'] = [f'name = ext.dat\nprint_every = {every}\ncol_1 = {{\ntype = force_energy\n}}']
            runs.append((side, os.path.join(base, side), ov))
    return runs


def compare_energy0(case, stage_dir, args):
    res = {}
    e, gpu = {}, {}
    run_dirs = sorted(d for d in os.listdir(stage_dir) if os.path.isdir(os.path.join(stage_dir, d)))
    for d in run_dirs:
        en, g = oxio.read_energy(os.path.join(stage_dir, d, 'energy.dat'), with_gpu=True)
        if len(en) == 0:
            return ERROR, f'{d}: no energy written', res
        e[d] = en[0, 1]
        if len(g) > 1 and np.isfinite(g[1]):
            gpu[d] = (en[1, 1], g[1])  # host and GPU energy of the same configuration
    stored = side_plan(case, 'old', args) is None
    if stored:
        st = stored_value(case, 'energy0')
        if st is None:
            return ERROR, 'no stored reference (run the CPU suite with --update-stored)', {'energies': e}
        for k, v in st.items():
            e['old' + k[3:]] = v
    res['energies'] = e
    if not all(np.isfinite(v) for v in e.values()):
        return FAIL, 'non-finite energy: ' + str(e), res
    tol_abs = case.get('energy0_tol', case['sys'].get('energy0_tol', 2e-6))
    cross_backend = side_plan(case, 'new', args)[1] != (side_plan(case, 'old', args) or [0, 'CPU'])[1]
    if cross_backend:
        tol_abs = max(tol_abs, 2e-5 + 1e-4 * abs(e['old']))
    msgs = []
    status = PASS
    for d in sorted(e):
        if d.startswith('new'):
            ref = 'old' + d[3:]
            diff = e[d] - e[ref]
            res[f'diff{d[3:]}'] = diff
            msgs.append(f'{d[3:] or "default"}: ref {e[ref]:.6f} test {e[d]:.6f} diff {diff:+.2e}'
                        + (' (stored ref)' if stored else ''))
            if abs(diff) > tol_abs:
                status = FAIL
    for d, (host, g) in sorted(gpu.items()):
        diff = g - host
        res[f'gpu_minus_host_{d}'] = diff
        tol_gpu = 2e-5 + 1e-4 * abs(host)
        msgs.append(f'{d}: GPU {g:.6f} vs host {host:.6f} ({diff:+.1e})')
        if abs(diff) > tol_gpu:
            status = FAIL
            msgs.append('GPU energy mismatch')
    # independent evaluation of the protein terms
    ref = case.get('reference')
    if ref and 'old' in e:
        conf = os.path.join(case['dir'], case['sys'].get('energy0_conf'))
        terms = anm_reference.energies(conf, os.path.join(case['dir'], ref['par']),
                                       case['sys'].get('protein_range'), ref.get('kb'), ref.get('kt'))
        N = terms.pop('N')
        prot = sum(terms.values()) / N
        res['reference_terms_per_particle'] = {k: v / N for k, v in terms.items()}
        res['reference_protein_energy'] = prot
        if case['sys'].get('nucleic_range') is None:
            for side in ('old', 'new'):
                d = e[side] - prot
                res[f'{side}_minus_reference'] = d
                msgs.append(f'{side}-numpy {d:+.2e}')
                if side == 'new' and abs(d) > max(tol_abs, 2e-6):
                    status = FAIL
        else:
            for side in ('old', 'new'):
                res[f'{side}_non_protein_terms'] = e[side] - prot
            msgs.append('non-protein terms (E - E_numpy_protein): ref %.6f test %.6f'
                        % (res['old_non_protein_terms'], res['new_non_protein_terms']))
    # Debye-Hueckel isolation for RNA hybrids
    if 'dh_off' in case.get('energy0_variants', {}) and 'old_dh_off' in e:
        dh_old = e['old'] - e['old_dh_off']
        dh_new = e['new'] - e['new_dh_off']
        res['DH_old'], res['DH_new'] = dh_old, dh_new
        ratio = dh_old / dh_new if dh_new != 0 else float('nan')
        res['DH_ratio_old_over_new'] = ratio
        msgs.append(f'DH term ref {dh_old:.6f} test {dh_new:.6f} ratio {ratio:.3f}')
    return status, '; '.join(msgs), res


def _min_image(d, box):
    return d - box * np.rint(d / box)


def compare_det(case, stage_dir, args):
    s = stage_settings(case, 'det', args)
    tol = s.get('tol', 1e-6)
    to, bo, xo = oxio.read_trajectory(os.path.join(stage_dir, 'old', 'trajectory.dat'), cols=15)
    tn, bn, xn = oxio.read_trajectory(os.path.join(stage_dir, 'new', 'trajectory.dat'), cols=15)
    if len(xo) == 0 or len(xn) == 0:
        return ERROR, 'empty trajectory', {}
    n = min(len(xo), len(xn))
    devs, vdevs = [], []
    for t in range(n):
        d = _min_image(xn[t][:, :3] - xo[t][:, :3], bo[t])
        devs.append(float(np.abs(d).max()))
        vdevs.append(float(np.abs(xn[t][:, 9:15] - xo[t][:, 9:15]).max()))
    eo = oxio.read_energy(os.path.join(stage_dir, 'old', 'energy.dat'))
    en = oxio.read_energy(os.path.join(stage_dir, 'new', 'energy.dat'))
    m = min(len(eo), len(en))
    ediff = float(np.abs(eo[:m, 1:3] - en[:m, 1:3]).max()) if m else float('nan')
    res = {'frames': n, 'max_pos_dev_per_frame': devs, 'max_vel_dev_per_frame': vdevs,
           'max_energy_diff': ediff, 'steps': s['steps']}
    first_bad = next((i for i, d in enumerate(devs) if d > tol), None)
    msg = (f'{n} frames over {s["steps"]} steps: max |dx| {max(devs):.1e} (final {devs[-1]:.1e}), '
           f'max |dv,dL| {max(vdevs):.1e}, max |dE| {ediff:.1e}')
    if tol == 0 and (max(vdevs) > 0 or ediff > 0):
        first_bad = first_bad if first_bad is not None else 0
    if first_bad is not None:
        msg += f'; exceeds {tol:g} from frame {first_bad} (t = {tn[first_bad]:g})'
        return FAIL, msg, res
    if len(xo) != len(xn):
        return FAIL, msg + f'; frame count differs ({len(xo)} vs {len(xn)})', res
    return PASS, msg + (' (bit-identical output)' if tol == 0 else ''), res


def _nve_one(rundir, ext_factor):
    en = oxio.read_energy(os.path.join(rundir, 'energy.dat'))
    if len(en) < 3 or not np.all(np.isfinite(en[:, 1:4])):
        raise NonFinite(f'{os.path.basename(rundir)}: non-finite or missing energies')
    t, U, K = en[:, 0], en[:, 1], en[:, 2]
    E = U + K
    ext_file = os.path.join(rundir, 'ext.dat')
    if os.path.isfile(ext_file):
        ext = oxio.read_columns(ext_file)[:, 0]
        m = min(len(ext), len(E))
        E, t, ext = E[:m], t[:m], ext[:m]
        E = E + ext_factor * ext
    slope, icpt = np.polyfit(t, E, 1)
    resid = E - (slope * t + icpt)
    return {'drift_per_time': float(slope), 'drift_total': float(slope * (t[-1] - t[0])),
            'fluct': float(resid.std()), 'E0': float(E[0]), 'max_dev': float(np.abs(E - E[0]).max()),
            'time': float(t[-1] - t[0])}


def compare_nve(case, stage_dir, args):
    s = stage_settings(case, 'nve', args)
    ext_factor = case['sys'].get('ext_energy_factor', 0.5)  # mutual traps come in pairs: force_energy double counts
    try:
        n = _nve_one(os.path.join(stage_dir, 'new'), ext_factor)
        if side_plan(case, 'old', args) is None:
            o = stored_value(case, 'nve')
            if o is None:
                return ERROR, 'no stored reference (run the CPU suite with --update-stored)', {'new': n}
        else:
            o = _nve_one(os.path.join(stage_dir, 'old'), ext_factor)
    except NonFinite as e:
        return FAIL, str(e), {}
    factor = s.get('factor', 3.0)
    d_floor = s.get('drift_floor', 2e-4)   # per particle, over the whole run
    f_floor = s.get('fluct_floor', 2e-4)
    res = {'old': o, 'new': n}
    msg = (f'conserved E drift over {n["time"]:g} time units: ref {o["drift_total"]:+.2e} test {n["drift_total"]:+.2e}; '
           f'fluct ref {o["fluct"]:.1e} test {n["fluct"]:.1e} (per particle)')
    status = PASS
    if abs(n['drift_total']) > factor * abs(o['drift_total']) + d_floor:
        status = FAIL
        msg += '; drift too large'
    if n['fluct'] > factor * o['fluct'] + f_floor:
        status = FAIL
        msg += '; fluctuations too large'
    return status, msg, res


class NonFinite(Exception):
    pass


def _stat_one(case, rundir, s):
    en = oxio.read_energy(os.path.join(rundir, 'energy.dat'))
    eq = s.get('equil_frac', 0.2)
    e = en[int(len(en) * eq):, 1]
    if len(e) == 0 or not np.all(np.isfinite(en[:, 1])):
        raise NonFinite(f'{os.path.basename(rundir)}: non-finite or missing energies')
    mean, se = analysis.block_average(e, s.get('energy_blocks', 10))
    out = {'E_mean': mean, 'E_se': se, 'n_energy': len(e)}
    t, b, x = oxio.read_trajectory(os.path.join(rundir, 'trajectory.dat'))
    if not np.all(np.isfinite(x)):
        raise NonFinite(f'{os.path.basename(rundir)}: non-finite coordinates in the trajectory')
    pr = case['sys'].get('protein_range') or [0, x.shape[1]]
    idx = np.arange(pr[0], pr[1])
    r, rse = analysis.trajectory_rmsf(x, b, idx, eq, s.get('rmsf_blocks', 5))
    out.update({'rmsf': r.tolist(), 'rmsf_se': rse.tolist(), 'n_frames': len(x)})
    nr = case['sys'].get('nucleic_range')
    if nr:
        r2, rse2 = analysis.trajectory_rmsf(x, b, np.arange(nr[0], nr[1]), eq, s.get('rmsf_blocks', 5))
        out.update({'rmsf_nucleic': r2.tolist(), 'rmsf_nucleic_se': rse2.tolist()})
    return out


def _stat_side(case, side_dir, s):
    """Average the replicas of one side. With several replicas the standard errors come from the
    spread of the independent replica means (robust to slow processes such as duplex fraying); with
    a single replica they come from block averaging."""
    reps = sorted(d for d in os.listdir(side_dir) if d.startswith('rep'))
    per = [_stat_one(case, os.path.join(side_dir, d), s) for d in reps]
    if len(per) == 1:
        return per[0]
    n = len(per)
    out = {'replicas': n, 'n_frames': sum(p['n_frames'] for p in per)}
    E = np.array([p['E_mean'] for p in per])
    out['E_mean'], out['E_se'] = float(E.mean()), float(E.std(ddof=1) / np.sqrt(n))
    out['E_replicas'] = E.tolist()
    for key in ('rmsf', 'rmsf_nucleic'):
        if key in per[0]:
            R = np.array([p[key] for p in per])
            out[key] = R.mean(0).tolist()
            out[key + '_se'] = (R.std(0, ddof=1) / np.sqrt(n)).tolist()
    return out


def compare_stat(case, stage_dir, args, variant=''):
    s = stage_settings(case, 'stat', args)
    suffix = f'_{variant}' if variant else ''
    try:
        n = _stat_side(case, os.path.join(stage_dir, 'new' + suffix), s)
        if side_plan(case, 'old', args) is None:
            o = stored_value(case, 'stat' + suffix)
            if o is None:
                return ERROR, 'no stored reference (run the CPU suite with --update-stored)', {'new': n}
        else:
            o = _stat_side(case, os.path.join(stage_dir, 'old' + suffix), s)
    except NonFinite as e:
        return FAIL, str(e), {}
    z_lim = s.get('z', 4.0)
    res = {'old': o, 'new': n}
    msgs = []
    status = PASS
    # potential energy
    dE = n['E_mean'] - o['E_mean']
    sig = np.hypot(o['E_se'], n['E_se'])
    floor = s.get('energy_floor', 2e-5)
    zE = dE / max(sig, 1e-12)
    res['E_diff'], res['E_sigma'], res['E_z'] = dE, sig, zE
    msgs.append(f'<E> old {o["E_mean"]:.5f}+-{o["E_se"]:.5f} new {n["E_mean"]:.5f}+-{n["E_se"]:.5f} (z {zE:+.1f})')
    if abs(dE) > s.get('z_energy', z_lim) * sig + floor:
        status = FAIL
        msgs.append('energy mismatch')
    # RMSF
    for key, label in (('rmsf', 'RMSF'), ('rmsf_nucleic', 'RMSF(NA)')):
        if key not in o:
            continue
        ro, rn = np.array(o[key]), np.array(n[key])
        so, sn = np.array(o[key + '_se']), np.array(n[key + '_se'])
        sig = np.hypot(so, sn)
        z = (rn - ro) / np.maximum(sig, 1e-6)
        n_out = int(np.sum(np.abs(z) > z_lim))
        frac = n_out / len(z)
        chi2 = float(np.mean(z ** 2))
        res[key + '_mean_z2'] = chi2
        mo, mn = ro.mean(), rn.mean()
        msig = np.hypot(so.mean(), sn.mean())  # conservative (fully correlated errors)
        res[key + '_frac_outliers'] = frac
        res[key + '_mean'] = (mo, mn)
        informational = key == 'rmsf_nucleic' and not s.get('test_nucleic_rmsf', False)
        msgs.append(f'{label} mean old {mo:.3f} new {mn:.3f} A, {n_out}/{len(z)} residues |z|>{z_lim:g}, '
                    f'<z^2> {chi2:.2f}'
                    + (' (not tested)' if informational else ''))
        if informational:
            continue
        # few replicas -> noisy standard errors -> heavy-tailed z: allow a few outliers, bound <z^2> and the mean
        if (n_out > max(1, s.get('max_outlier_frac', 0.05) * len(z)) or chi2 > s.get('max_mean_z2', 4.0)
                or abs(mn - mo) > z_lim * msig + 0.02 * mo):
            status = FAIL
            msgs.append(f'{label} mismatch')
    # ANM prediction (informational)
    msf_file = os.path.join(case['dir'], 'anm_msf.dat')
    if os.path.isfile(msf_file) and case['sys'].get('nucleic_range') is None:
        msf = np.loadtxt(msf_file)
        pred = np.sqrt(msf[:, 3])
        rn = np.array(n['rmsf'])
        if len(pred) == len(rn):
            corr = float(np.corrcoef(pred, rn)[0, 1])
            res['anm_prediction'] = {'pearson_new': corr, 'ratio_new_over_anm': float(rn.mean() / pred.mean()),
                                     'pearson_old': float(np.corrcoef(pred, np.array(o['rmsf']))[0, 1])}
            msgs.append(f'vs ANM prediction: r={corr:.2f}, sim/ANM={rn.mean() / pred.mean():.2f}')
    return status, '; '.join(msgs), res


def check_inputs(case):
    """The original code takes the protein bonds from the topology neighbour lists (plus the
    C-side neighbour) while the current code takes them from the parameter file: make sure they agree."""
    inp = dict(case['sys'].get('common', {}))
    inp.update(case.get('common', {}))
    inp.update(case.get('old', {}))
    top, par = inp.get('topology'), inp.get('parfile')
    if top is None or par is None or str(par) == 'none':
        return SKIP, 'no topology/parfile', {}
    with open(os.path.join(case['dir'], top)) as f:
        if '5->3' in f.readline():
            return SKIP, 'new-style topology (bonds come from the par file only)', {}
    header, parts = oxio.read_old_topology(os.path.join(case['dir'], top))
    springs = set(tuple(sorted(s)) for s in oxio.read_par(os.path.join(case['dir'], par)))
    mixed = len(header) >= 5
    bonds = set()
    for i, (strand, code, n3, n5, extra) in enumerate(parts):
        if mixed and strand > 0:
            continue
        for j in extra + ([n5] if n5 >= 0 else []):
            bonds.add(tuple(sorted((i, j))))
    missing = springs - bonds
    extra = bonds - springs
    # extra topology bonds without a spring are harmless (k = 0) as long as they are backbone links
    msg = f'{len(springs)} springs, {len(bonds)} topology bonds, {len(missing)} springs not in topology'
    if missing:
        return FAIL, msg + f' (e.g. {sorted(missing)[:3]})', {'extra': len(extra)}
    return PASS, msg, {'extra': len(extra)}



# ----------------------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--old', required=True, help='reference binary: the original anm-oxdna (default mode) or a '
                                                 'CPU build of the code under test (--same-code)')
    ap.add_argument('--new', required=True, help='oxDNA binary under test')
    ap.add_argument('--backend', default='CPU', choices=['CPU', 'CUDA'], help='backend of the --new runs')
    ap.add_argument('--old-backend', default='CPU', choices=['CPU', 'CUDA'], help='backend of the --old runs')
    ap.add_argument('--cuda-precision', default='mixed', choices=['mixed', 'float', 'double'],
                    help='backend_precision of the CUDA runs (double needs a -DCUDA_DOUBLE=ON build)')
    ap.add_argument('--same-code', action='store_true',
                    help='both sides use the "new" inputs of each case (CPU vs CUDA of the same oxDNA, or oxNA)')
    ap.add_argument('--old-inputs', default=None, choices=['old', 'new'], help='deprecated alias: new = --same-code')
    ap.add_argument('--quick', action='store_true', help='only the quick subset, with shorter runs')
    ap.add_argument('--smoke', action='store_true', help='only the smoke subset (energy0 + a short nve), fail fast')
    ap.add_argument('--only', nargs='*', help='glob(s) on "system/case" or "system"')
    ap.add_argument('--stages', default=','.join(ALL_STAGES), help='comma separated subset of ' + ','.join(ALL_STAGES))
    ap.add_argument('--jobs', '-j', type=int, default=os.cpu_count(), help='simultaneous simulations')
    ap.add_argument('--gpu-jobs', type=int, default=4, help='simultaneous CUDA simulations (they share the GPU)')
    ap.add_argument('--out', default=None, help='output directory (default results/<backend>[_quick])')
    ap.add_argument('--timeout', type=float, default=3600)
    ap.add_argument('--update-stored', action='store_true',
                    help='write the test-side results of "stored" cases to systems/<s>/stored.json')
    ap.add_argument('--list', action='store_true', help='list the selected cases and exit')
    args = ap.parse_args()
    if args.old_inputs == 'new':
        args.same_code = True
    args.mode = 'same_code' if args.same_code else 'original'
    tag = lambda b: b if b == 'CPU' else ('CUDA' if args.cuda_precision == 'mixed' else f'CUDA_{args.cuda_precision}')
    args.btag, args.old_btag = tag(args.backend), tag(args.old_backend)
    args.old, args.new = os.path.abspath(args.old), os.path.abspath(args.new)
    if args.out is None:
        sub = 'smoke' if args.smoke else ('quick' if args.quick else '')
        args.out = os.path.join(ROOT, 'results', args.btag + ('_same' if args.same_code else '') + (f'_{sub}' if sub else ''))
    args.out = os.path.abspath(args.out)
    args.stages = args.stages.split(',')
    if args.smoke:
        args.stages = [s for s in args.stages if s in ('energy0', 'nve')]
    return args


def det_possible(case, args):
    o = side_plan(case, 'old', args)
    return o is not None and o[1] == side_plan(case, 'new', args)[1]


def main():
    args = parse_args()
    cases = load_cases(args)
    if args.list:
        for c in cases:
            print(c['name'], '-', c.get('description', ''))
        return 0
    if not cases:
        print('no cases selected')
        return 1
    if os.path.isdir(args.out):
        shutil.rmtree(args.out)
    os.makedirs(args.out)

    # 1. all simulations, in parallel (CUDA runs additionally limited by --gpu-jobs)
    jobs = []
    for c in cases:
        for st in enabled_stages(c, args):
            if st == 'det' and not det_possible(c, args):
                continue
            for side, rundir, ov in stage_runs(c, st, args):
                binary, btag, sets = side_plan(c, side, args)
                make_input(c, sets, btag, ov, rundir)
                cost = max(ov.get('steps', 0), 1) * c['sys'].get('cost', 1.0)
                jobs.append((cost, binary, rundir, c['name'], st, side, btag.startswith('CUDA')))
    jobs.sort(key=lambda j: -j[0])  # longest first
    print(f'{len(cases)} cases, {len(jobs)} simulations, {args.jobs} parallel jobs -> {args.out}', flush=True)
    t0 = time.time()
    runinfo = {}
    gpu_sem = threading.Semaphore(args.gpu_jobs)
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_oxdna, j[1], j[2], args.timeout, gpu_sem if j[6] else None): j for j in jobs}
        for k, fut in enumerate(cf.as_completed(futs)):
            j = futs[fut]
            r = fut.result()
            runinfo[j[2]] = r
            tag = 'ok' if r['error'] is None else 'ERROR: ' + r['error']
            print(f'[{k + 1}/{len(jobs)} {time.time() - t0:6.0f}s] {j[3]} {j[4]} {os.path.relpath(j[2], os.path.join(args.out, j[3], j[4]))} '
                  f'({r["time"]:.0f}s) {tag}', flush=True)

    # 2. analysis
    results = []
    for c in cases:
        cres = {'case': c['name'], 'description': c.get('description', ''), 'checks': {}}
        if args.mode == 'original' and c.get('old_code', 'original') == 'original':
            st, msg, data = check_inputs(c)
            cres['checks']['inputs'] = {'status': st, 'message': msg, 'data': data}
        for stg in enabled_stages(c, args):
            if stg == 'det' and not det_possible(c, args):
                why = 'not comparable across backends' if side_plan(c, 'old', args) else 'stored reference'
                cres['checks'][stg] = {'status': SKIP, 'message': why}
                continue
            stage_dir = os.path.join(args.out, c['name'], stg)
            errs = [f'{os.path.relpath(d, stage_dir)}: {r["error"]}' for d, r in runinfo.items()
                    if d.startswith(stage_dir + os.sep) and r['error']]
            try:
                if errs:
                    status, msg, data = ERROR, '; '.join(sorted(errs)[:4]) + (' ...' if len(errs) > 4 else ''), {}
                elif stg == 'energy0':
                    status, msg, data = compare_energy0(c, stage_dir, args)
                elif stg == 'det':
                    status, msg, data = compare_det(c, stage_dir, args)
                elif stg == 'nve':
                    status, msg, data = compare_nve(c, stage_dir, args)
                else:
                    status, msg, data = compare_stat(c, stage_dir, args)
                    for vname in c.get('stat_variants', {}):
                        vs, vm, vd = compare_stat(c, stage_dir, args, vname)
                        cres['checks'][f'stat_{vname}'] = {'status': vs, 'message': vm, 'data': vd}
            except Exception as e:  # noqa: BLE001
                status, msg, data = ERROR, f'analysis failed: {e!r}', {}
            cres['checks'][stg] = {'status': status, 'message': msg, 'data': data}
        if args.update_stored and c.get('old_code') == 'stored' and args.mode == 'original':
            _store(c, cres, args)
        # expectations
        for name, chk in cres['checks'].items():
            if name == 'inputs' or chk['status'] == SKIP:
                continue
            x = expectation(c, 'stat' if name.startswith('stat_') else name, args)
            if x is None:
                continue
            applies, reason = _xfail_applies(x, chk)
            chk['expected_failure'] = reason
            if reason:
                if chk['status'] in (FAIL, ERROR) and not applies:
                    chk['expected_failure'] = None
                elif chk['status'] in (FAIL, ERROR):
                    chk['status'] = XFAIL
                elif chk['status'] == PASS:
                    chk['status'] = XPASS
        results.append(cres)

    # 3. summary
    mode = 'same code (reference = CPU build of the code under test)' if args.same_code else 'original vs current'
    lines = [f'# ANM/oxNA regression results ({args.btag}{" quick" if args.quick else ""}{" smoke" if args.smoke else ""})',
             '', f'- mode: {mode}',
             f'- old/reference: `{args.old}` ({args.old_btag})', f'- new/test: `{args.new}` ({args.btag})',
             f'- wall time: {time.time() - t0:.0f} s', '',
             '| case | check | status | details |', '|---|---|---|---|']
    counts = {}
    for r in results:
        for name, chk in r['checks'].items():
            det = chk['message'].replace('|', '/')
            if chk.get('expected_failure'):
                det += f' [expected: {chk["expected_failure"]}]'
            lines.append(f'| {r["case"]} | {name} | {chk["status"]} | {det} |')
            counts[chk['status']] = counts.get(chk['status'], 0) + 1
    nfail = counts.get(FAIL, 0) + counts.get(ERROR, 0)
    lines += ['', 'totals: ' + ', '.join(f'{k} {v}' for k, v in sorted(counts.items())), '',
              f'**{nfail} unexpected failure(s)**']
    summary = '\n'.join(lines) + '\n'
    with open(os.path.join(args.out, 'summary.md'), 'w') as f:
        f.write(summary)
    with open(os.path.join(args.out, 'results.json'), 'w') as f:
        json.dump(results, f, indent=1, default=float)
    print()
    print(summary)
    return 1 if nfail else 0


def _store(case, cres, args):
    info = {'binary': args.new, 'date': time.strftime('%Y-%m-%d')}
    chk = cres['checks']
    if 'energy0' in chk:
        e = chk['energy0'].get('data', {}).get('energies', {})
        vals = {k: v for k, v in e.items() if k.startswith('new')}
        if vals:
            update_stored(case, 'energy0', vals, info)
    for name in [n for n in chk if n.startswith('stat') or n == 'nve']:
        d = chk[name].get('data', {})
        if 'new' in d:
            new = d['new']
            if name.startswith('stat'):
                new = {k: new[k] for k in ('E_mean', 'E_se', 'rmsf', 'rmsf_se', 'rmsf_nucleic', 'rmsf_nucleic_se')
                       if k in new}
            update_stored(case, name, new, info)


if __name__ == '__main__':
    sys.exit(main())
