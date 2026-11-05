#!/usr/bin/env python3
"""Turn an ANM parameter file into an ANMT one by appending the four reference cosines to every
backbone spring (j == i + 1), computed from the orientations stored in a configuration:

    a0 = r_ij . a1_i,  b0 = -r_ij . a1_j,  c0 = a1_i . a1_j,  d0 = a3_i . a3_j   (r_ij unit vector, min. image)

This is the same definition used by ANMUtils/models.py, so the angular energy of CONF is zero.
Optionally appends per-spring kb kt. Used to build ANMT inputs for systems whose ANM network was
generated earlier (e.g. the hybrid cage).
"""
import argparse
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('conf')
    ap.add_argument('par_in')
    ap.add_argument('par_out')
    ap.add_argument('--kb', type=float)
    ap.add_argument('--kt', type=float)
    a = ap.parse_args()
    with open(a.conf) as f:
        lines = f.readlines()
    box = np.array([float(x) for x in lines[1].split()[2:5]])
    data = np.array([[float(x) for x in l.split()] for l in lines[3:] if l.strip()])
    pos, a1, a3 = data[:, :3], data[:, 3:6], data[:, 6:9]
    a1 = a1 / np.linalg.norm(a1, axis=1)[:, None]
    a3 = a3 - a1 * (a1 * a3).sum(1)[:, None]
    a3 = a3 / np.linalg.norm(a3, axis=1)[:, None]
    out = []
    with open(a.par_in) as f:
        out.append(f.readline())
        for line in f:
            s = line.split()
            if len(s) < 5:
                continue
            i, j = int(s[0]), int(s[1])
            if j - i == 1:
                r = pos[j] - pos[i]
                r -= box * np.rint(r / box)
                r /= np.linalg.norm(r)
                cos = [r @ a1[i], -r @ a1[j], a1[i] @ a1[j], a3[i] @ a3[j]]
                s = s[:5] + ['%.15g' % c for c in cos]
                if a.kb is not None:
                    s += ['%g' % a.kb, '%g' % a.kt]
            out.append(' '.join(s) + '\n')
    with open(a.par_out, 'w') as f:
        f.writelines(out)


if __name__ == '__main__':
    main()
