#!/usr/bin/env python3
"""Independent numpy evaluation of the protein (ANM / ANMT) energy terms.

This does not use either oxDNA code. It reads a configuration, a .par file and
(optionally) the protein index range, and returns per-term energies (totals, not
per particle):

  spring  : sum over par-file springs of 0.5 k (r - r0)^2
  excl    : protein-protein excluded volume over ALL protein pairs (bonded or not)
  angular : ANMT bending + torsion over backbone springs (j == i + 1) with angle params

Positions are used as written in the configuration with minimum-image distances.
"""
import argparse
import sys
import numpy as np

EXCL_EPS = 2.0
SIGMA, RSTAR, B, RC = 0.35, 0.349, 306484596.421, 0.352894


def read_conf(path):
    with open(path) as f:
        lines = f.readlines()
    box = np.array([float(x) for x in lines[1].split()[2:5]])
    data = np.array([[float(x) for x in l.split()] for l in lines[3:] if l.strip()])
    return box, data


def read_par(path):
    springs, angles = [], {}
    with open(path) as f:
        f.readline()
        for line in f:
            s = line.split()
            if len(s) < 5:
                continue
            i, j, r0, k = int(s[0]), int(s[1]), float(s[2]), float(s[4])
            springs.append((i, j, r0, k))
            if j - i == 1 and len(s) >= 9:
                kbkt = (float(s[9]), float(s[10])) if len(s) >= 11 else None
                angles[i] = (tuple(float(x) for x in s[5:9]), kbkt)
    return springs, angles


def excl(r):
    r = np.asarray(r)
    e = np.zeros_like(r)
    quart = (r < RC) & (r > RSTAR)
    lj = r <= RSTAR
    e[quart] = EXCL_EPS * B * (r[quart] - RC) ** 4
    l6 = (SIGMA / r[lj]) ** 6
    e[lj] = 4 * EXCL_EPS * (l6 ** 2 - l6)
    return e


def minimg(d, box):
    return d - box * np.rint(d / box)


def energies(conf, par, protein=None, kb=None, kt=None):
    box, data = read_conf(conf)
    pos, a1, a3 = data[:, 0:3], data[:, 3:6], data[:, 6:9]
    springs, angles = read_par(par)
    out = {"spring": 0.0, "excl": 0.0, "angular": 0.0}
    for i, j, r0, k in springs:
        r = np.linalg.norm(minimg(pos[j] - pos[i], box))
        out["spring"] += 0.5 * k * (r - r0) ** 2
    if protein is None:
        idx = np.arange(len(pos))
    else:
        idx = np.arange(protein[0], protein[1])
    P = pos[idx]
    for n in range(len(P) - 1):
        d = np.linalg.norm(minimg(P[n + 1:] - P[n], box), axis=1)
        out["excl"] += excl(d[d < RC]).sum()
    for i, ((a0, b0, c0, d0), kbkt) in angles.items():
        j = i + 1
        kbend, ktor = kbkt if kbkt is not None else (kb, kt)
        if kbend is None:
            raise ValueError("bending/torsion constants not given")
        rij = minimg(pos[j] - pos[i], box)
        rij /= np.linalg.norm(rij)
        o1 = rij @ a1[i] - a0
        o2 = -rij @ a1[j] - b0
        o3 = a1[i] @ a1[j] - c0
        o4 = a3[i] @ a3[j] - d0
        out["angular"] += 0.5 * kbend * (o1 ** 2 + o2 ** 2) + 0.5 * ktor * (o3 ** 2 + o4 ** 2)
    out["N"] = len(pos)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("conf")
    ap.add_argument("par")
    ap.add_argument("--protein", nargs=2, type=int, help="first, last+1 index of the protein block")
    ap.add_argument("--kb", type=float)
    ap.add_argument("--kt", type=float)
    a = ap.parse_args()
    e = energies(a.conf, a.par, a.protein, a.kb, a.kt)
    N = e.pop("N")
    for k, v in e.items():
        print(f"{k:8s} total {v:.10g}  per-particle {v / N:.10g}")
