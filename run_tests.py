#!/usr/bin/env python3
"""Regression tests for the oxDNA protein (ANM) models: original implementation vs current oxDNA.

Each case (systems/<system>/config.json -> "cases") is simulated with a reference binary (--old) and a
test binary (--new) and compared through up to three stages:

  energy0  potential energy of a perturbed configuration (steps = 0); also checked against an
           independent numpy evaluation of the protein terms (tools/anm_reference.py)
  det      short MD run with the same seed and thermostat: trajectories must be (near) identical
  stat     long MD runs with different seeds: block-averaged potential energy and per-residue RMSF
           must agree within statistical error

Example:
  ./run_tests.py --old OLD/bin/oxDNA --new NEW/bin/oxDNA --backend CPU --quick
"""
import argparse
import concurrent.futures as cf
import fnmatch
import json
import os
import shutil
import subprocess
import sys
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
BACKEND_INPUT = {
    'CPU': {'backend': 'CPU', 'backend_precision': 'double'},
    'CUDA': {'backend': 'CUDA', 'backend_precision': 'mixed', 'CUDA_list': 'verlet', 'CUDA_sort_every': 0,
             'use_edge': 1, 'edge_n_forces': 1, 'max_density_multiplier': 10},
}
FILE_KEYS = ('topology', 'parfile', 'conf_file', 'external_forces_file', 'massfile')

PASS, FAIL, XFAIL, XPASS, SKIP, ERROR = 'PASS', 'FAIL', 'XFAIL', 'XPASS', 'SKIP', 'ERROR'


# ----------------------------------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------------------------------

def load_cases(only, quick):
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
            if only and not any(fnmatch.fnmatch(full, pat) or fnmatch.fnmatch(name, pat) for pat in only):
                continue
            if quick and not case.get('quick', False):
                continue
            c = dict(case)
            c['system'] = name
            c['name'] = full
            c['dir'] = os.path.join(sysdir, name)
            c['sys'] = cfg
            cases.append(c)
    return cases


def stage_settings(case, stage, quick):
    s = dict(case['sys'].get('stages', {}).get(stage, {}))
    s.update(case.get('stages', {}).get(stage, {}))
    if quick:
        s.update(s.pop('quick', {}))
    else:
        s.pop('quick', None)
    return s


def expectation(case, stage, backend):
    """Return the reason string if this stage is expected to fail on this backend, else None."""
    for key in (f'{stage}@{backend}', f'*@{backend}', stage, '*'):
        if key in case.get('xfail', {}):
            return case['xfail'][key]
    return None


def make_input(case, side, backend, stage_overrides, rundir):
    inp = dict(BASE_INPUT)
    inp.update(BACKEND_INPUT[backend])
    inp.update(case['sys'].get('common', {}))
    inp.update(case.get('common', {}))
    inp.update(case.get(side, {}))
    inp.update(case.get(f'{side}@{backend}', {}))
    inp.update(stage_overrides)
    for k in FILE_KEYS:
        if k in inp and not os.path.isabs(str(inp[k])):
            inp[k] = os.path.join(case['dir'], inp[k])
    inp.update({'trajectory_file': 'trajectory.dat', 'lastconf_file': 'last_conf.dat',
                'energy_file': 'energy.dat', 'log_file': 'log.dat'})
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(rundir, 'input'), 'w') as f:
        for k, v in inp.items():
            f.write(f'{k} = {v}\n')
    return inp


def run_oxdna(binary, rundir, timeout):
    t0 = time.time()
    with open(os.path.join(rundir, 'stdout.txt'), 'w') as out:
        try:
            p = subprocess.run([binary, 'input'], cwd=rundir, stdout=out, stderr=subprocess.STDOUT,
                               timeout=timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    elapsed = time.time() - t0
    with open(os.path.join(rundir, 'stdout.txt')) as f:
        txt = f.read()
    err = None
    if rc != 0 or 'ERROR' in txt:
        lines = [l for l in txt.splitlines() if 'ERROR' in l]
        err = lines[-1].strip() if lines else f'exit code {rc}'
    return {'rc': rc, 'time': elapsed, 'error': err}


# ----------------------------------------------------------------------------------------------------
# stages
# ----------------------------------------------------------------------------------------------------

def stage_runs(case, stage, args):
    """Return the list of (side, rundir, input-overrides) for a stage."""
    s = stage_settings(case, stage, args.quick)
    base = os.path.join(args.out, case['name'], stage)
    runs = []
    if stage == 'energy0':
        ov = {'steps': 0, 'print_energy_every': 1, 'print_conf_interval': 1000, 'seed': 1,
              'conf_file': case['sys'].get('energy0_conf', case['sys']['common']['conf_file'])}
        for side in ('old', 'new'):
            runs.append((side, os.path.join(base, side), ov))
            for vname, vov in case.get('energy0_variants', {}).items():
                o = dict(ov)
                o.update(vov)
                runs.append((side, os.path.join(base, f'{side}_{vname}'), o))
    elif stage == 'det':
        n = s['steps']
        ov = {'steps': n, 'print_energy_every': s.get('print_every', n // 5),
              'print_conf_interval': s.get('print_every', n // 5), 'seed': 1}
        for side in ('old', 'new'):
            runs.append((side, os.path.join(base, side), ov))
    elif stage == 'stat':
        for side, seed in (('old', 1000), ('new', 2000)):
            for vname, vov in [('', {})] + list(case.get('stat_variants', {}).items()):
                for rep in range(s.get('replicas', 1)):
                    ov = {'steps': s['steps'], 'print_energy_every': s['print_energy_every'],
                          'print_conf_interval': s['print_conf_interval'], 'seed': seed + rep}
                    if 'stat_conf' in case['sys']:
                        ov['conf_file'] = case['sys']['stat_conf']
                    ov.update(vov)
                    runs.append((side, os.path.join(base, side + (f'_{vname}' if vname else ''), f'rep{rep}'), ov))
    return runs


def compare_energy0(case, stage_dir, args):
    res = {}
    e = {}
    for d in sorted(os.listdir(stage_dir)):
        en = oxio.read_energy(os.path.join(stage_dir, d, 'energy.dat'))
        if len(en) == 0:
            return ERROR, f'{d}: no energy written', res
        e[d] = en[0, 1]
    res['energies'] = e
    if not all(np.isfinite(v) for v in e.values()):
        return FAIL, 'non-finite energy: ' + str(e), res
    tol_abs = case['sys'].get('energy0_tol', 2e-6)
    if args.backend != args.old_backend:
        tol_abs = max(tol_abs, 1e-4 * max(abs(e['old']), 1.0))
    msgs = []
    status = PASS
    for d in e:
        if d.startswith('new'):
            ref = 'old' + d[3:]
            diff = e[d] - e[ref]
            res[f'diff{d[3:]}'] = diff
            msgs.append(f'{d[3:] or "default"}: old {e[ref]:.6f} new {e[d]:.6f} diff {diff:+.2e}')
            if d == 'new' and abs(diff) > tol_abs:
                status = FAIL
    # independent evaluation of the protein terms
    ref = case.get('reference')
    if ref:
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
                msgs.append(f'{side}-ref {d:+.2e}')
                if side == 'new' and abs(d) > tol_abs:
                    status = FAIL
        else:
            for side in ('old', 'new'):
                res[f'{side}_non_protein_terms'] = e[side] - prot
            msgs.append('non-protein terms (E - E_ref_protein): old %.6f new %.6f'
                        % (res['old_non_protein_terms'], res['new_non_protein_terms']))
    # optional variant analysis (e.g. the Debye-Hueckel isolation for RNA hybrids)
    if 'dh_off' in case.get('energy0_variants', {}):
        dh_old = e['old'] - e['old_dh_off']
        dh_new = e['new'] - e['new_dh_off']
        res['DH_old'], res['DH_new'] = dh_old, dh_new
        ratio = dh_old / dh_new if dh_new != 0 else float('nan')
        res['DH_ratio_old_over_new'] = ratio
        msgs.append(f'DH term old {dh_old:.6f} new {dh_new:.6f} ratio {ratio:.3f}')
        diff_off = e['new_dh_off'] - e['old_dh_off']
        if abs(diff_off) > tol_abs:
            msgs.append('MISMATCH with DH switched off')
            status = FAIL
    return status, '; '.join(msgs), res


def _min_image(d, box):
    return d - box * np.rint(d / box)


def compare_det(case, stage_dir, args):
    s = stage_settings(case, 'det', args.quick)
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
    if first_bad is not None:
        msg += f'; exceeds {tol:g} from frame {first_bad} (t = {tn[first_bad]:g})'
        return FAIL, msg, res
    if len(xo) != len(xn):
        return FAIL, msg + f'; frame count differs ({len(xo)} vs {len(xn)})', res
    return PASS, msg, res


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
    s = stage_settings(case, 'stat', args.quick)
    suffix = f'_{variant}' if variant else ''
    try:
        o = _stat_side(case, os.path.join(stage_dir, 'old' + suffix), s)
        n = _stat_side(case, os.path.join(stage_dir, 'new' + suffix), s)
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
    if top is None or par is None:
        return SKIP, 'no topology/parfile', {}
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

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--old', required=True, help='reference oxDNA binary (original ANM implementation)')
    ap.add_argument('--new', required=True, help='oxDNA binary under test')
    ap.add_argument('--backend', default='CPU', choices=['CPU', 'CUDA'], help='backend of the --new runs')
    ap.add_argument('--old-backend', default='CPU', choices=['CPU', 'CUDA'], help='backend of the --old runs')
    ap.add_argument('--old-inputs', default='old', choices=['old', 'new'],
                    help="which input set of each case the --old runs use: 'new' compares two builds of the "
                         "same code (e.g. CPU vs CUDA of the current oxDNA)")
    ap.add_argument('--quick', action='store_true', help='only the quick subset, with shorter runs')
    ap.add_argument('--only', nargs='*', help='glob(s) on "system/case" or "system"')
    ap.add_argument('--stages', default='energy0,det,stat', help='comma separated subset of energy0,det,stat')
    ap.add_argument('--jobs', '-j', type=int, default=os.cpu_count())
    ap.add_argument('--out', default=None, help='output directory (default results/<backend>[_quick])')
    ap.add_argument('--timeout', type=float, default=3600)
    ap.add_argument('--list', action='store_true', help='list the selected cases and exit')
    args = ap.parse_args()
    args.old, args.new = os.path.abspath(args.old), os.path.abspath(args.new)
    if args.out is None:
        args.out = os.path.join(ROOT, 'results', args.backend + ('_quick' if args.quick else ''))
    args.out = os.path.abspath(args.out)
    stages = args.stages.split(',')

    cases = load_cases(args.only, args.quick)
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

    # 1. all simulations, in parallel
    jobs = []
    for c in cases:
        for st in stages:
            if st not in c.get('stages_enabled', ['energy0', 'det', 'stat']):
                continue
            if st == 'det' and (args.backend != args.old_backend):
                continue  # trajectories cannot be compared across backends
            for side, rundir, ov in stage_runs(c, st, args):
                backend = args.old_backend if side == 'old' else args.backend
                binary = args.old if side == 'old' else args.new
                make_input(c, args.old_inputs if side == 'old' else 'new', backend, ov, rundir)
                cost = ov.get('steps', 0) * c['sys'].get('cost', 1.0)
                jobs.append((cost, binary, rundir, c['name'], st, side))
    jobs.sort(key=lambda j: -j[0])  # longest first
    print(f'{len(cases)} cases, {len(jobs)} simulations, {args.jobs} parallel jobs -> {args.out}', flush=True)
    t0 = time.time()
    runinfo = {}
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_oxdna, j[1], j[2], args.timeout): j for j in jobs}
        for k, fut in enumerate(cf.as_completed(futs)):
            j = futs[fut]
            r = fut.result()
            runinfo[j[2]] = r
            tag = 'ok' if r['error'] is None else 'ERROR: ' + r['error']
            print(f'[{k + 1}/{len(jobs)} {time.time() - t0:6.0f}s] {j[3]} {j[4]} '
                  f'{os.path.basename(j[2])} ({r["time"]:.0f}s) {tag}', flush=True)

    # 2. analysis
    results = []
    for c in cases:
        cres = {'case': c['name'], 'description': c.get('description', ''), 'checks': {}}
        st, msg, data = check_inputs(c)
        cres['checks']['inputs'] = {'status': st, 'message': msg, 'data': data}
        for stg in stages:
            if stg not in c.get('stages_enabled', ['energy0', 'det', 'stat']):
                continue
            if stg == 'det' and args.backend != args.old_backend:
                cres['checks'][stg] = {'status': SKIP, 'message': 'not comparable across backends'}
                continue
            stage_dir = os.path.join(args.out, c['name'], stg)
            errs = [f'{os.path.basename(d)}: {r["error"]}' for d, r in runinfo.items()
                    if d.startswith(stage_dir + os.sep) and r['error']]
            try:
                if errs:
                    status, msg, data = ERROR, '; '.join(errs), {}
                elif stg == 'energy0':
                    status, msg, data = compare_energy0(c, stage_dir, args)
                elif stg == 'det':
                    status, msg, data = compare_det(c, stage_dir, args)
                else:
                    status, msg, data = compare_stat(c, stage_dir, args)
                    for vname in c.get('stat_variants', {}):
                        vs, vm, vd = compare_stat(c, stage_dir, args, vname)
                        cres['checks'][f'stat_{vname}'] = {'status': vs, 'message': vm, 'data': vd}
            except Exception as e:  # noqa: BLE001
                status, msg, data = ERROR, f'analysis failed: {e!r}', {}
            cres['checks'][stg] = {'status': status, 'message': msg, 'data': data}
        # expectations
        for name, chk in cres['checks'].items():
            if name == 'inputs':
                continue
            reason = expectation(c, name.split('_')[0] if name.startswith('stat_') else name, args.backend)
            chk['expected_failure'] = reason
            if reason:
                if chk['status'] in (FAIL, ERROR):
                    chk['status'] = XFAIL
                elif chk['status'] == PASS:
                    chk['status'] = XPASS
        results.append(cres)

    # 3. summary
    lines = [f'# ANM regression results ({args.backend}{" quick" if args.quick else ""})', '',
             f'- old: `{args.old}` ({args.old_backend})', f'- new: `{args.new}` ({args.backend})',
             f'- wall time: {time.time() - t0:.0f} s', '',
             '| case | check | status | details |', '|---|---|---|---|']
    nfail = 0
    for r in results:
        for name, chk in r['checks'].items():
            det = chk['message'].replace('|', '/')
            if chk.get('expected_failure'):
                det += f' [expected: {chk["expected_failure"]}]'
            lines.append(f'| {r["case"]} | {name} | {chk["status"]} | {det} |')
            if chk['status'] in (FAIL, ERROR):
                nfail += 1
    lines += ['', f'**{nfail} unexpected failure(s)**']
    summary = '\n'.join(lines) + '\n'
    with open(os.path.join(args.out, 'summary.md'), 'w') as f:
        f.write(summary)
    with open(os.path.join(args.out, 'results.json'), 'w') as f:
        json.dump(results, f, indent=1, default=float)
    print()
    print(summary)
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
