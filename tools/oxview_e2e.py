#!/usr/bin/env python3
"""End-to-end check of oxView's protein/ANM exports against oxDNA and oat.

Input: a folder written by oxView (scratch/oxview_e2e): one sub-folder per case with system.{top,dat,par}
(new-style topology), old_topology/system.{top,dat,par}, oxview_info.json and, for fitted cases,
fluctuation_fit.json; the source PDBs are in inputs/.

Checks per case:
  files    every .par spring has r0 = |r_j - r_i| (minimum image) in the exported configuration, every amino-acid
           pair with C-alphas within the cutoff has a spring, ANMT cosines match the exported frames
  oat      oat anm_fluctuation_solver run on the source PDB with oxView's settings gives the same network
           (pairs, r0, k, cosines), configuration (positions, frames) and, for fitted cases, the same k and msf
  energy   oxDNA reads both topologies for each intended interaction; steps=0 energies of the new- and
           old-style exports agree; on a perturbed configuration the protein energy of protein-only systems
           agrees with an independent numpy evaluation (tools/anm_reference.py); wrapping every particle into
           the box (so the protein crosses the boundary) leaves the energy unchanged
  md       short thermostatted MD (finite energies) followed by a constant-energy run (drift per particle)

  oxview_e2e.py <oxview folder> <oxDNA binary> [--out DIR] [--only CASE ...]
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anm_reference as ref  # noqa: E402

A_TO_SIM = 1 / 8.518

CASES = {
    '01_1ubq_anm': ('1UBQ.pdb', 'anm', ['AC', 'DNANM']),
    '02_1ubq_anmt': ('1UBQ.pdb', 'anmt', ['ACT', 'DNACT']),
    '03a_1aay_protein_dna_anm': ('1AAY.pdb', 'anm', ['DNANM']),
    '03b_1aay_protein_dna_anmt': ('1AAY.pdb', 'anmt', ['DNACT']),
    '04a_1urn_protein_rna_anm': ('1URN.pdb', 'anm', ['RNANM']),
    '04b_1urn_protein_rna_anmt': ('1URN.pdb', 'anmt', ['RNACT']),
    '05_kdpg_box_anm': ('kdpg.pdb', 'anm', ['AC', 'DNANM']),
    '06a_1ubq_fit_anm': ('1UBQ.pdb', 'anm', ['AC']),
    '06b_1ubq_fit_anmt': ('1UBQ.pdb', 'anmt', ['ACT']),
}

BASE_INPUT = {
    'sim_type': 'MD', 'backend': 'CPU', 'backend_precision': 'double', 'T': '300K', 'verlet_skin': 0.2,
    'salt_concentration': 1.0, 'dt': 0.002, 'thermostat': 'john', 'newtonian_steps': 51, 'diff_coeff': 1,
    'time_scale': 'linear', 'restart_step_counter': 1, 'no_stdout_energy': 1, 'refresh_vel': 1,
    'energy_file': 'energy.dat', 'log_file': 'log.dat', 'trajectory_file': 'trajectory.dat',
    'lastconf_file': 'last_conf.dat',
}


class Result:
    def __init__(self):
        self.rows = []

    def add(self, case, check, ok, msg):
        status = 'PASS' if ok else 'FAIL'
        self.rows.append((case, check, status, msg))
        print(f'{status} {case} {check}: {msg}', flush=True)


def protein_mask(top_path):
    """Boolean mask of the amino acids, from either topology style."""
    with open(top_path) as f:
        header = f.readline().split()
        lines = [l for l in f if l.strip()]
    if header[-1] == '5->3':
        mask = []
        for l in lines:
            seq, *kw = l.split()
            kw = dict(k.split('=', 1) for k in kw if '=' in k)
            n = len(seq)  # one-letter sequences only in these exports
            mask += [kw.get('type') == 'peptide'] * n
        return np.array(mask)
    return np.array([int(l.split()[0]) < 0 for l in lines])


def check_files(case, d, cutoff_A, res):
    box, data = ref.read_conf(os.path.join(d, 'system.dat'))
    springs, angles = ref.read_par(os.path.join(d, 'system.par'))
    prot = np.where(protein_mask(os.path.join(d, 'system.top')))[0]
    pos, a1, a3 = data[:, 0:3], data[:, 3:6], data[:, 6:9]
    with open(os.path.join(d, 'system.par')) as f:
        n_header = int(f.readline().split()[0])
    ok = n_header == len(prot)
    msgs = [f'{len(springs)} springs, {len(prot)} amino acids (par header {n_header})']
    ij = np.array([(s[0], s[1]) for s in springs])
    r0 = np.array([s[2] for s in springs])
    r = np.linalg.norm(ref.minimg(pos[ij[:, 1]] - pos[ij[:, 0]], box), axis=1)
    dr = np.abs(r - r0).max()
    msgs.append(f'max |r - r0| {dr:.1e}')
    ok &= dr < 1e-5
    if not np.all(np.isin(ij, prot)):
        ok = False
        msgs.append('springs on non-protein particles')
    # every protein pair within the cutoff must have a spring, none beyond it
    P = pos[prot]
    dd = np.linalg.norm(ref.minimg(P[:, None] - P[None], box), axis=2)
    iu = np.triu_indices(len(P), 1)
    within = {(prot[a], prot[b]) for a, b, x in zip(*iu, dd[iu]) if x <= cutoff_A * A_TO_SIM}
    have = {tuple(sorted(p)) for p in ij.tolist()}
    missing, extra = within - have, have - within
    msgs.append(f'pairs within {cutoff_A} A: {len(within)}, missing {len(missing)}, beyond cutoff {len(extra)}')
    ok &= not missing and not extra
    if angles:
        worst = 0.
        for i, ((c0, c1, c2, c3), kbkt) in angles.items():
            j = i + 1
            u = ref.minimg(pos[j] - pos[i], box)
            u /= np.linalg.norm(u)
            calc = (u @ a1[i], -u @ a1[j], a1[i] @ a1[j], a3[i] @ a3[j])
            worst = max(worst, max(abs(x - y) for x, y in zip(calc, (c0, c1, c2, c3))))
        n_bb = sum(1 for i, j in ij.tolist() if j == i + 1)
        msgs.append(f'ANMT: {len(angles)} of {n_bb} backbone springs with cosines, max cosine error {worst:.1e}')
        ok &= worst < 1e-5 and len(angles) == n_bb
    res.add(case, 'files', ok, '; '.join(msgs))
    return prot


def check_oat(case, d, pdb, model, cutoff_A, prot, res):
    from oxDNA_analysis_tools import anm_fluctuation_solver as afs
    st = afs.read_structure(pdb)
    fit_file = os.path.join(d, 'fluctuation_fit.json')
    info = json.load(open(os.path.join(d, 'oxview_info.json')))
    nw = info['network']
    kw = dict(model=model, cutoff=float(cutoff_A), temperature=300.)
    if model == 'anmt':
        kw.update(kb=float(nw.get('kb', 50)), kt=float(nw.get('kt', 50)), write_kbkt=True)
    fitted = os.path.isfile(fit_file)
    if fitted:
        # oxView: one global scale, no-intercept slope 1 over all residues = oat's linear fit without outliers
        kw.update(max_deviations=1e12)
    else:
        kw.update(spring_constant=float(nw.get('k_unfitted_sim', 1)))
    net = afs.anm_fluctuation_solver(st, **kw)
    msgs, ok = [], True
    if st.n_residues != len(prot):
        res.add(case, 'oat', False, f'oat reads {st.n_residues} residues, oxView exported {len(prot)}')
        return
    # oxView's particle index of every oat residue (protein particles in order)
    springs, angles = ref.read_par(os.path.join(d, 'system.par'))
    ov = {(min(i, j), max(i, j)): (r0, k) for i, j, r0, k in springs}
    oat = {(min(prot[i], prot[j]), max(prot[i], prot[j])): (r0, k) for i, j, r0, k in
           zip(net.i.tolist(), net.j.tolist(), net.r0.tolist(), net.k.tolist())}
    common = set(ov) & set(oat)
    msgs.append(f'springs oxView {len(ov)} oat {len(oat)} common {len(common)}')
    ok &= len(common) == len(ov) == len(oat)
    if common:
        dr0 = max(abs(ov[p][0] - oat[p][0]) for p in common)
        dk = max(abs(ov[p][1] / oat[p][1] - 1) for p in common)
        msgs.append(f'max |dr0| {dr0:.1e}, max |k_ov/k_oat - 1| {dk:.1e} (k oat {np.mean([v[1] for v in oat.values()]):.4f})')
        ok &= dr0 < 1e-5 and dk < (1e-3 if fitted else 1e-9)
    # configuration: positions up to a translation, frames
    box, data = ref.read_conf(os.path.join(d, 'system.dat'))
    ca = st.coords * A_TO_SIM
    x = data[prot, 0:3]
    shift = (x - ca).mean(axis=0)
    dpos = np.abs(x - ca - shift).max()
    a1e, a3e = afs._engine_frames(net.a1s, net.a3s)
    da1 = np.abs(data[prot, 3:6] - a1e).max()
    da3 = np.abs(data[prot, 6:9] - a3e).max()
    msgs.append(f'positions (after centring) max diff {dpos:.1e}, a1 {da1:.1e}, a3 {da3:.1e}')
    ok &= dpos < 1e-5 and da1 < 1e-4 and da3 < 1e-4
    if model == 'anmt' and net.angles is not None:
        bb = np.where(net.j == net.i + 1)[0]
        worst = 0.
        for p, a in zip(bb, net.angles):
            i = prot[net.i[p]]
            if i in angles:
                worst = max(worst, np.abs(np.array(angles[i][0]) - a).max())
        msgs.append(f'cosines max diff {worst:.1e}')
        ok &= worst < 1e-4
    if fitted:
        fit = json.load(open(fit_file))
        msf_oat = net.pred_bfactors / afs.BFACTOR_CONVERSION
        msf_ov = np.array(fit['fitted_msf_A2'])
        rel = np.abs(msf_ov / msf_oat - 1).max()
        msgs.append(f'mean k oxView {fit["mean_k_sim"]:.4f} oat {net.k.mean():.4f}; fitted msf max rel diff {rel:.1e}')
        ok &= abs(fit['mean_k_sim'] / net.k.mean() - 1) < 1e-3 and rel < 1e-3
    res.add(case, 'oat', ok, '; '.join(msgs))


def write_input(path, extra):
    with open(path, 'w') as f:
        for k, v in {**BASE_INPUT, **extra}.items():
            f.write(f'{k} = {v}\n')


def run_oxdna(binary, rundir, extra):
    os.makedirs(rundir, exist_ok=True)
    write_input(os.path.join(rundir, 'input'), extra)
    p = subprocess.run([binary, 'input'], cwd=rundir, capture_output=True, text=True, timeout=3600)
    if p.returncode != 0:
        log = os.path.join(rundir, 'log.dat')
        text = p.stdout + p.stderr + (open(log).read() if os.path.isfile(log) else '')
        err = [l for l in text.splitlines() if 'ERROR' in l]
        raise RuntimeError(err[-1] if err else f'exit code {p.returncode}')
    return np.loadtxt(os.path.join(rundir, 'energy.dat'), ndmin=2)


def write_conf(path, box, data):
    with open(path, 'w') as f:
        f.write('t = 0\n')
        f.write('b = ' + ' '.join(f'{x:.17g}' for x in box) + '\n')
        f.write('E = 0 0 0\n')
        for row in data:
            f.write(' '.join(f'{x:.17g}' for x in row) + '\n')


def check_energy_md(case, d, interactions, prot, out, binary, res):
    nested = {'new': d, 'old': os.path.join(d, 'old_topology')}
    n_part = len(ref.read_conf(os.path.join(d, 'system.dat'))[1])
    # nucleic acids converted from crystal structures can sit just outside the FENE range (and PDB gaps give
    # stretched bonds): cap the backbone force, as for any relaxation of a converted structure
    mbf = {'max_backbone_force': 5, 'max_backbone_force_far': 10} if len(prot) < n_part else {}
    for inter in interactions:
        e0 = {}
        try:
            for style, sd in nested.items():
                files = {'topology': os.path.join(sd, 'system.top'), 'parfile': os.path.join(sd, 'system.par'),
                         'interaction_type': inter, **mbf}
                rd = os.path.join(out, case, inter, style)
                e0[style] = run_oxdna(binary, os.path.join(rd, 'energy0'), {
                    **files, 'conf_file': os.path.join(sd, 'system.dat'), 'steps': 0, 'print_energy_every': 1,
                    'print_conf_interval': 1000, 'seed': 1})[0, 1]
                # every particle wrapped into the box: the protein (and the NA) cross the boundary
                box, data = ref.read_conf(os.path.join(sd, 'system.dat'))
                wrapped = data.copy()
                wrapped[:, 0:3] = np.mod(data[:, 0:3] + 0.37 * box, box)  # shifted so that the cut is off-centre
                write_conf(os.path.join(rd, 'wrapped.dat'), box, wrapped)
                ew = run_oxdna(binary, os.path.join(rd, 'energy0_wrapped'), {
                    **files, 'conf_file': os.path.join(rd, 'wrapped.dat'), 'steps': 0, 'print_energy_every': 1,
                    'print_conf_interval': 1000, 'seed': 1})[0, 1]
                crossing = np.any(np.abs(ref.minimg(wrapped[:, None, :3] - wrapped[None, :, :3], box) -
                                         (wrapped[:, None, :3] - wrapped[None, :, :3])) > 1e-9)
                res.add(case, f'{inter} wrap ({style})', abs(ew - e0[style]) < 1e-6,
                        f'E unwrapped {e0[style]:.8f} wrapped {ew:.8f} (particles split across the box: {bool(crossing)})')
            res.add(case, f'{inter} new vs old top', abs(e0['new'] - e0['old']) < 1e-6,
                    f'E new {e0["new"]:.8f} old {e0["old"]:.8f}')

            # perturbed configuration of the new-style export: protein energy vs numpy (protein-only systems)
            box, data = ref.read_conf(os.path.join(d, 'system.dat'))
            rng = np.random.default_rng(7)
            pert = data.copy()
            pert[:, 0:3] += rng.normal(0, 0.02, (len(data), 3))
            rd = os.path.join(out, case, inter, 'new')
            write_conf(os.path.join(rd, 'perturbed.dat'), box, pert)
            files = {'topology': os.path.join(d, 'system.top'), 'parfile': os.path.join(d, 'system.par'),
                     'interaction_type': inter, **mbf}
            ep = run_oxdna(binary, os.path.join(rd, 'energy0_perturbed'), {
                **files, 'conf_file': os.path.join(rd, 'perturbed.dat'), 'steps': 0, 'print_energy_every': 1,
                'print_conf_interval': 1000, 'seed': 1})[0, 1]
            if len(prot) == len(data):
                r = ref.energies(os.path.join(rd, 'perturbed.dat'), os.path.join(d, 'system.par'), (0, len(data)))
                en = (r['spring'] + r['excl'] + r['angular']) / r['N']
                res.add(case, f'{inter} perturbed vs numpy', abs(ep - en) < 1e-6 * max(1, abs(en)),
                        f'oxDNA {ep:.8f} numpy {en:.8f} (spring {r["spring"] / r["N"]:.5f} '
                        f'excl {r["excl"] / r["N"]:.5f} angular {r["angular"] / r["N"]:.5f})')

            # MD: 20 time units thermostatted, then 10 time units at constant energy
            md = run_oxdna(binary, os.path.join(rd, 'md'), {
                **files, 'conf_file': os.path.join(d, 'system.dat'), 'dt': 0.002, 'steps': 10000,
                'print_energy_every': 500, 'print_conf_interval': 10000, 'seed': 11})
            nve = run_oxdna(binary, os.path.join(rd, 'nve'), {
                **files, 'conf_file': os.path.join(rd, 'md', 'last_conf.dat'), 'dt': 0.001, 'steps': 10000,
                'thermostat': 'no', 'refresh_vel': 0, 'print_energy_every': 100, 'print_conf_interval': 10000,
                'seed': 12})
            finite = np.all(np.isfinite(md)) and np.all(np.isfinite(nve))
            E = nve[:, 1] + nve[:, 2]
            slope = np.polyfit(nve[:, 0], E, 1)[0] * (nve[-1, 0] - nve[0, 0])
            fluct = np.std(E - np.polyval(np.polyfit(nve[:, 0], E, 1), nve[:, 0]))
            Tk = md[len(md) // 2:, 2].mean()
            res.add(case, f'{inter} md+nve', finite and abs(slope) < 2e-3 and fluct < 2e-3,
                    f'MD U {md[-1, 1]:+.4f} K {Tk:.4f}; NVE drift {slope:+.1e} fluct {fluct:.1e} per particle over 10 time units')
        except Exception as ex:  # noqa: BLE001
            res.add(case, f'{inter}', False, f'oxDNA error: {ex}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder')
    ap.add_argument('oxdna')
    ap.add_argument('--out', default=None)
    ap.add_argument('--only', nargs='*')
    a = ap.parse_args()
    folder = os.path.abspath(a.folder)
    out = os.path.abspath(a.out or os.path.join(folder, 'check'))
    binary = os.path.abspath(a.oxdna)
    res = Result()
    for case, (pdb, model, inters) in CASES.items():
        if a.only and case not in a.only:
            continue
        d = os.path.join(folder, case)
        if not os.path.isdir(d):
            res.add(case, 'present', False, 'missing')
            continue
        cutoff = json.load(open(os.path.join(d, 'oxview_info.json')))['network']['cutoff_A']
        prot = check_files(case, d, cutoff, res)
        try:
            check_oat(case, d, os.path.join(folder, 'inputs', pdb), model, cutoff, prot, res)
        except Exception as ex:  # noqa: BLE001
            res.add(case, 'oat', False, f'error: {ex}')
        check_energy_md(case, d, inters, prot, out, binary, res)
    lines = ['| case | check | status | details |', '|---|---|---|---|']
    lines += [f'| {c} | {k} | {s} | {m} |' for c, k, s, m in res.rows]
    n_fail = sum(r[2] != 'PASS' for r in res.rows)
    lines.append(f'\n{len(res.rows) - n_fail} passed, {n_fail} failed')
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'summary.md'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines[-1:]))
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
