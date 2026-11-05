#!/usr/bin/env python3
"""Generate an ideal A-form oxRNA duplex (old-style, 3'->5' topology).

Python 3 port of StrandGenerator.a_generate from the oxDNA UTILS/generate-RNA.py script
(rise 0.3287, twist 32.7 deg, inclination 15.5 deg, diameter 2.35, CM_CENTER_DS = RNA_POS_BASE = 0.4).

Usage: make_rna_duplex.py SEQUENCE OUT_PREFIX [--box 20]
writes OUT_PREFIX.dat and OUT_PREFIX.top
"""
import argparse
import math
import numpy as np

COMP = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


def rotz(deg):
    t = math.radians(deg)
    return np.array([[math.cos(t), -math.sin(t), 0], [math.sin(t), math.cos(t), 0], [0, 0, 1.]])


def a_generate(n, inclination=15.5, diameter=2.35, rise=0.3287, rot=32.7, cm_center_ds=0.4):
    bb = 2.0
    cord = math.cos(math.radians(inclination)) * bb
    c2c = math.sqrt((diameter / 2.) ** 2 - (cord / 2.) ** 2)
    r1 = np.array([c2c, +cord / 2., -(bb / 2.) * math.sin(math.radians(inclination))])
    r2 = np.array([c2c, -cord / 2., +(bb / 2.) * math.sin(math.radians(inclination))])
    Z = np.array([0., 0., 1.])
    R = rotz(rot)
    s1, s2 = [], []
    ci, si = math.cos(math.radians(inclination)), math.sin(math.radians(inclination))
    for _ in range(n):
        r1 = R @ r1 + rise * Z
        r2 = R @ r2 + rise * Z
        a1 = (r2 - r1) / np.linalg.norm(r2 - r1)
        a1p = np.array([a1[0], a1[1], 0.])
        a3 = -(-ci * Z + si * a1p / np.linalg.norm(a1p))
        s1.append((r1 + cm_center_ds * a1, a1, a3))
    RT = R.T
    for _ in range(n):
        a1 = -(r2 - r1) / np.linalg.norm(r2 - r1)
        a1p = np.array([a1[0], a1[1], 0.])
        a3 = -(ci * Z + si * a1p / np.linalg.norm(a1p))
        s2.append((r2 + cm_center_ds * a1, a1, a3))
        r1 = RT @ r1 - rise * Z
        r2 = RT @ r2 - rise * Z
    return s1, s2


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('sequence', help='sequence of strand 1 (the second strand is its complement)')
    ap.add_argument('out')
    ap.add_argument('--box', type=float, default=20.0)
    a = ap.parse_args()
    seq1 = a.sequence.upper().replace('T', 'U')
    seq2 = ''.join(COMP[c] for c in reversed(seq1))
    s1, s2 = a_generate(len(seq1))
    n = len(seq1)
    with open(a.out + '.dat', 'w') as f:
        f.write('t = 0\nb = {0} {0} {0}\nE = 0 0 0\n'.format(a.box))
        for pos, a1, a3 in s1 + s2:
            vals = list(pos) + list(a1) + list(a3) + [0.] * 6
            f.write(' '.join('%.15g' % v for v in vals) + '\n')
    with open(a.out + '.top', 'w') as f:
        f.write(f'{2 * n} 2\n')
        for s, seq in enumerate((seq1, seq2)):
            for i, c in enumerate(seq):
                idx = s * n + i
                n3 = idx - 1 if i > 0 else -1
                n5 = idx + 1 if i < n - 1 else -1
                f.write(f'{s + 1} {c} {n3} {n5}\n')


if __name__ == '__main__':
    main()
