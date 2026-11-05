#!/usr/bin/env python3
"""Combine a protein system (from make_protein_system.py) and a nucleic-acid duplex into a hybrid system.

The protein comes first in the topology (required by the original DNACT/RNACT parameter-file checks,
which reject particle ids larger than the number of amino acids). The nucleic acid is translated so that
its closest nucleotide sits --gap simulation units from the closest amino acid, and a pair of mutual
traps (written in the old '{}' external-forces format, read by both codes) keeps the two in contact.

Writes OUTDIR/conf.dat, hybrid.top, anm.par, anmt.par and forces.txt.
"""
import argparse
import os
import shutil
import numpy as np


def read_conf(path):
    with open(path) as f:
        lines = f.readlines()
    box = float(lines[1].split()[2])
    data = np.array([[float(x) for x in l.split()] for l in lines[3:] if l.strip()])
    if data.shape[1] < 15:
        data = np.hstack([data, np.zeros((len(data), 15 - data.shape[1]))])
    return box, data


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

    # put the nucleic acid far along d, then slide it back until the closest pair is at distance gap
    pcom, ncom = prot[:, :3].mean(0), nuc[:, :3].mean(0)
    nuc[:, :3] += pcom - ncom + 30 * d
    def mind(n):
        diff = n[:, None, :3] - prot[None, :, :3]
        return np.linalg.norm(diff, axis=2).min()
    shift = 0.0
    while mind(nuc) > a.gap:
        nuc[:, :3] -= 0.01 * d
    diff = nuc[:, None, :3] - prot[None, :, :3]
    dist = np.linalg.norm(diff, axis=2)
    inuc, ipro = np.unravel_index(dist.argmin(), dist.shape)
    r0 = dist[inuc, ipro]

    npro = len(prot)
    allc = np.vstack([prot, nuc])
    allc[:, :3] -= allc[:, :3].mean(0) - a.box / 2
    with open(os.path.join(a.outdir, 'conf.dat'), 'w') as f:
        f.write('t = 0\nb = {0} {0} {0}\nE = 0 0 0\n'.format(a.box))
        for row in allc:
            f.write(' '.join('%.15g' % v for v in row) + '\n')

    with open(os.path.join(a.protein_dir, 'anm.top')) as f:
        plines = f.read().splitlines()
    nps = int(plines[0].split()[1])
    with open(a.nucleic_top) as f:
        nlines = f.read().splitlines()
    nnuc, nns = (int(x) for x in nlines[0].split()[:2])
    with open(os.path.join(a.outdir, 'hybrid.top'), 'w') as f:
        f.write(f'{npro + nnuc} {nps + nns} {nnuc} {npro} {nns}\n')
        for l in plines[1:]:
            f.write(l + '\n')
        for l in nlines[1:]:
            s, b, n3, n5 = l.split()[:4]
            n3, n5 = int(n3), int(n5)
            n3 = n3 + npro if n3 >= 0 else -1
            n5 = n5 + npro if n5 >= 0 else -1
            f.write(f'{s} {b} {n3} {n5}\n')
    for par in ('anm.par', 'anmt.par'):
        shutil.copy(os.path.join(a.protein_dir, par), os.path.join(a.outdir, par))
    with open(os.path.join(a.outdir, 'forces.txt'), 'w') as f:
        for p, q in ((ipro, npro + inuc), (npro + inuc, ipro)):
            f.write('{\ntype = mutual_trap\nparticle = %d\nref_particle = %d\nstiff = %g\nr0 = %.4f\nPBC = 1\n}\n\n'
                    % (p, q, a.trap_stiff, r0))
    print(f'protein {npro}, nucleotides {nnuc}, closest pair {ipro}-{npro + inuc} at {r0:.3f}')


if __name__ == '__main__':
    main()
