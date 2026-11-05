#!/usr/bin/env python3
"""Generate ANM and ANMT simulation files for a PDB with the (patched) original ANMUtils models.py.

Usage:
    python tools/make_protein_system.py PDB OUTDIR [--cutoff 13] [--kb KB --kt KT]

Writes into OUTDIR:
    conf.dat          initial configuration (CA positions, a1 = CA->side-chain COM, a3 perpendicular)
    anm.par           ANM springs (i j r0 s k)
    anmt.par          ANMT springs: backbone (j == i+1) lines carry 4 cosines (+ kb kt if --kb/--kt given)
    anm.top           protein-only old-style topology ('N N_strands', with neighbour lists), for AC/ACT
    mixed.top         same topology with the mixed header 'N N_strands 0 N 0', for DNANM/DNACT/RNANM/RNACT
    anm_msf.dat       residue index, experimental B-factor, ANM-predicted B-factor, ANM-predicted MSF (A^2)

Both ANM and ANMT use the same unitary-ANM fit to the B-factors (as in the original notebooks: ANMT
B-factors cannot be computed analytically, so ANMT reuses the fitted ANM spring constant).
"""
import argparse
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models as m  # noqa: E402  (patched copy, see tools/README_models.md)


def export(model, pdb, workdir):
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        m.export_to_simulation(model, pdb, upstreamdir='')
    finally:
        os.chdir(cwd)
    return [os.path.join(workdir, 'generated.' + ext) for ext in ('par', 'top', 'dat')]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pdb')
    ap.add_argument('outdir')
    ap.add_argument('--cutoff', type=float, default=13.0, help='ANM cutoff in Angstrom (default 13)')
    ap.add_argument('--kb', type=float, default=0.0, help='per-spring bending constant written to anmt.par')
    ap.add_argument('--kt', type=float, default=0.0, help='per-spring torsion constant written to anmt.par')
    ap.add_argument('--min-box', type=float, default=20.0, help='minimum box side (simulation units)')
    args = ap.parse_args()

    pdb = os.path.abspath(args.pdb)
    os.makedirs(args.outdir, exist_ok=True)
    coord, bfact, bmap = m.get_pdb_info(pdb, returntype='cbm')

    anm = m.ANM(coord, bfact, T=300, cutoff=args.cutoff)
    anm.calc_ANM_unitary(cuda=False)
    anmt = m.ANMT(coord, bfact, T=300, cutoff=args.cutoff, kb=args.kb, kt=args.kt)
    anmt.calc_ANM_unitary(cuda=False)

    with tempfile.TemporaryDirectory() as tmp:
        par, top, dat = export(anm, pdb, tmp)
        shutil.copy(par, os.path.join(args.outdir, 'anm.par'))
        shutil.copy(top, os.path.join(args.outdir, 'anm.top'))
        shutil.copy(dat, os.path.join(args.outdir, 'conf.dat'))
        par, top, dat = export(anmt, pdb, tmp)
        shutil.copy(par, os.path.join(args.outdir, 'anmt.par'))
        # ANMT and ANM share the same configuration and network topology
        for new, ref in ((dat, 'conf.dat'), (top, 'anm.top')):
            with open(new) as a, open(os.path.join(args.outdir, ref)) as b:
                if a.read() != b.read():
                    raise RuntimeError(f'ANM and ANMT exports produced different {ref}')

    # models.py sets the box to 2.5 * span, which is too small for short peptides (the neighbour lists of
    # the hybrid interactions need a box larger than twice the DNA/RNA cutoff)
    conf = os.path.join(args.outdir, 'conf.dat')
    with open(conf) as f:
        clines = f.readlines()
    box = float(clines[1].split()[2])
    if box < args.min_box:
        clines[1] = 'b = {0} {0} {0}\n'.format(args.min_box)
        with open(conf, 'w') as f:
            f.writelines(clines)

    # mixed-header topology for the hybrid interactions (protein only: 0 nucleotides, 0 nucleic strands)
    with open(os.path.join(args.outdir, 'anm.top')) as f:
        lines = f.read().splitlines()
    N, nstr = lines[0].split()[:2]
    with open(os.path.join(args.outdir, 'mixed.top'), 'w') as f:
        f.write(f'{N} {nstr} 0 {N} 0\n')
        f.write('\n'.join(lines[1:]) + '\n')

    with open(os.path.join(args.outdir, 'anm_msf.dat'), 'w') as f:
        f.write('# index exp_bfactor anm_bfactor anm_msf_A2   (ANM gamma = %.6g pN/A, k_sim = %.6g)\n'
                % (anm.ana_gamma, anm.ana_gamma / 0.05709))
        for i, (be, ba, msd) in enumerate(zip(anm.exp_bfactors, anm.ana_bfactors, anm.ana_msd)):
            f.write(f'{i} {be:.6g} {ba:.6g} {msd:.6g}\n')
    print('wrote', args.outdir)


if __name__ == '__main__':
    main()
