#!/usr/bin/env python3
"""Make systems/ubq_dna_wall/wall.dat: the ubq_dna configuration with amino acids inside the protein-DNA walls.

The protein-nucleotide excluded volume acts between an amino acid and the backbone and base sites of a nucleotide. It
is a repulsive Lennard-Jones potential that turns into a quartic below r_c, scaled by a stiffness: 1 for both walls in
the CPU code, while master's CUDA DNANM used EXCL_EPS = 2 (oxDNA PR #240). In equilibrated configurations amino acids
are rarely inside a wall, so the regression energies did not test it. This moves the amino acid closest to a backbone
site into the quartic part of that wall (r* < r < r_c) and the amino acid closest to a base site into the LJ part
(r < r*), and prints every protein-DNA pair inside a wall with its energy (stiffness 1).

  make_ubq_dna_wall.py            (writes systems/ubq_dna_wall/wall.dat)
"""
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, 'systems', 'ubq_dna', 'equilibrated.dat')
DST = os.path.join(REPO, 'systems', 'ubq_dna_wall', 'wall.dat')
PROTEIN, NUCLEIC = range(0, 76), range(76, 100)

# oxDNA2 (DNANM) sites relative to the nucleotide centre, and the protein-DNA walls (ANMModel.cpp)
POS_MM_BACK1, POS_MM_BACK2, POS_BASE = -0.34, 0.3408, 0.4
WALLS = {
    'backbone': dict(sigma=0.57, rstar=0.569, b=178699253.5, rc=0.572934, stiffness=1.0, target=0.5695),
    'base': dict(sigma=0.36, rstar=0.359, b=296866090.0, rc=0.362897, stiffness=1.0, target=0.35),
}


def wall_energy(w, r):
    if r >= w['rc']:
        return 0.
    if r > w['rstar']:
        return w['stiffness'] * w['b'] * (r - w['rc']) ** 4
    lj = (w['sigma'] / r) ** 6
    return 4 * w['stiffness'] * (lj ** 2 - lj)


def read(path):
    lines = open(path).read().splitlines()
    box = np.array([float(x) for x in lines[1].split()[2:5]])
    data = np.array([[float(x) for x in l.split()] for l in lines[3:] if l.strip()])
    return lines[:3], box, data


def sites(data):
    pos, a1, a3 = data[:, 0:3], data[:, 3:6], data[:, 6:9]
    a2 = np.cross(a3, a1)
    return {'backbone': pos + POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2, 'base': pos + POS_BASE * a1}


def minimg(d, box):
    return d - box * np.rint(d / box)


def contacts(data, box):
    s = sites(data)
    out = []
    for kind, w in WALLS.items():
        for j in NUCLEIC:
            d = np.linalg.norm(minimg(data[list(PROTEIN), 0:3] - s[kind][j], box), axis=1)
            for i in np.where(d < w['rc'])[0]:
                out.append((kind, int(i), j, float(d[i]), wall_energy(w, d[i])))
    return out


def main():
    head, box, data = read(SRC)
    moved = set()
    for kind, w in WALLS.items():
        s = sites(data)[kind]
        best = None
        for j in NUCLEIC:
            d = minimg(data[list(PROTEIN), 0:3] - s[j], box)
            r = np.linalg.norm(d, axis=1)
            for i in np.argsort(r):
                if int(i) in moved:
                    continue
                if best is None or r[i] < best[0]:
                    best = (r[i], int(i), j, d[i])
                break
        r, i, j, d = best
        shift = d / r * (w['target'] - r)
        data[i, 0:3] += shift
        moved.add(i)
        print(f'{kind}: amino acid {i} moved by {np.linalg.norm(shift):.3f} to {w["target"]} from the {kind} site of nucleotide {j} '
              f'(r* = {w["rstar"]}, r_c = {w["rc"]})')

    total = 0.
    print('protein-DNA pairs inside a wall:')
    for kind, i, j, r, e in contacts(data, box):
        total += e
        print(f'  {kind:8s} amino acid {i:2d} - nucleotide {j:3d}: r = {r:.4f}, E = {e:.6f}')
    n = len(data)
    print(f'total wall energy {total:.6f} ({total / n:.2e} per particle): a wall twice as strong would change the energy per '
          f'particle by {total / n:.2e}')
    os.makedirs(os.path.dirname(DST), exist_ok=True)
    with open(DST, 'w') as f:
        f.write('\n'.join(head) + '\n')
        for row in data:
            f.write(' '.join(f'{x:.10g}' for x in row) + '\n')
    print(f'written to {DST}')


if __name__ == '__main__':
    main()
