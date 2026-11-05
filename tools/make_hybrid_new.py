#!/usr/bin/env python3
"""Like make_hybrid.py, but for a nucleic-acid system with a new-style (5'->3', type=DNA/RNA) topology.
The result has a new-style topology (protein strands first, type=peptide), as NANM/NACT require.

    make_hybrid_new.py PROTEIN_DIR NUCLEIC_CONF NUCLEIC_TOP OUTDIR [--gap 1.0] [--box 30]

Writes OUTDIR/conf.dat, topology.top, anm.par, anmt.par and forces.txt (a mutual-trap pair between
the closest amino acid and nucleotide).
"""
import argparse
import os
import shutil

import numpy as np

from make_hybrid import read_conf


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('protein_dir')
    ap.add_argument('nucleic_conf')
    ap.add_argument('nucleic_top')
    ap.add_argument('outdir')
    ap.add_argument('--gap', type=float, default=1.0)
    ap.add_argument('--box', type=float, default=30.0)
    ap.add_argument('--direction', nargs=3, type=float, default=[1, 0, 0])
    ap.add_argument('--trap-stiff', type=float, default=1.0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    _, prot = read_conf(os.path.join(a.protein_dir, 'conf.dat'))
    _, nuc = read_conf(a.nucleic_conf)
    d = np.array(a.direction) / np.linalg.norm(a.direction)
    nuc[:, :3] += prot[:, :3].mean(0) - nuc[:, :3].mean(0) + 30 * d
    dist = lambda: np.linalg.norm(nuc[:, None, :3] - prot[None, :, :3], axis=2)
    while dist().min() > a.gap:
        nuc[:, :3] -= 0.01 * d
    D = dist()
    inuc, ipro = np.unravel_index(D.argmin(), D.shape)
    r0 = D[inuc, ipro]
    npro = len(prot)
    allc = np.vstack([prot, nuc])
    allc[:, :3] -= allc[:, :3].mean(0) - a.box / 2
    with open(os.path.join(a.outdir, 'conf.dat'), 'w') as f:
        f.write('t = 0\nb = {0} {0} {0}\nE = 0 0 0\n'.format(a.box))
        for row in allc:
            f.write(' '.join('%.15g' % v for v in row) + '\n')
    # peptide strands from the old-style protein topology (N to C order)
    with open(os.path.join(a.protein_dir, 'anm.top')) as f:
        f.readline()
        chains = {}
        for l in f:
            s = l.split()
            if s:
                chains.setdefault(int(s[0]), []).append(s[1])
    with open(a.nucleic_top) as f:
        h = f.readline().split()
        nlines = [l.strip() for l in f if l.strip()]
    nnuc = int(h[0])
    with open(os.path.join(a.outdir, 'topology.top'), 'w') as f:
        f.write(f'{npro + nnuc} {len(chains) + len(nlines)} 5->3\n')
        for k in sorted(chains, key=abs):
            f.write(''.join(chains[k]) + ' type=peptide\n')
        for l in nlines:
            f.write(l + '\n')
    for par in ('anm.par', 'anmt.par'):
        shutil.copy(os.path.join(a.protein_dir, par), os.path.join(a.outdir, par))
    with open(os.path.join(a.outdir, 'forces.txt'), 'w') as f:
        for p, q in ((ipro, npro + inuc), (npro + inuc, ipro)):
            f.write('{\ntype = mutual_trap\nparticle = %d\nref_particle = %d\nstiff = %g\nr0 = %.4f\nPBC = 1\n}\n\n'
                    % (p, q, a.trap_stiff, r0))
    print(f'protein {npro}, nucleotides {nnuc}, closest pair {ipro}-{npro + inuc} at {r0:.3f}')


if __name__ == '__main__':
    main()
