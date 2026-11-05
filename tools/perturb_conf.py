#!/usr/bin/env python3
"""Write a randomly perturbed copy of an oxDNA configuration (deterministic for a given seed).

Positions get Gaussian noise of width --sigma; a1/a3 are rotated by a random rotation whose angle is
Gaussian with width --sigma-rot (radians). Used to make initial-energy tests non-trivial: the
configurations written by the ANM generator sit exactly at the network minimum.
"""
import argparse
import numpy as np


def rotation(axis, angle):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('conf')
    ap.add_argument('out')
    ap.add_argument('--sigma', type=float, default=0.02)
    ap.add_argument('--sigma-rot', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=12345)
    ap.add_argument('--range', nargs=2, type=int, help='only perturb particles first..last-1')
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    with open(a.conf) as f:
        lines = f.readlines()
    head, body = lines[:3], [l for l in lines[3:] if l.strip()]
    first, last = (a.range if a.range else (0, len(body)))
    out = list(head)
    for i, l in enumerate(body):
        v = np.array([float(x) for x in l.split()])
        if first <= i < last:
            v[0:3] += rng.normal(0, a.sigma, 3)
            R = rotation(rng.normal(size=3), rng.normal(0, a.sigma_rot))
            a1 = R @ v[3:6]
            a3 = R @ v[6:9]
            a3 -= a1 * (a1 @ a3) / (a1 @ a1)
            v[3:6] = a1 / np.linalg.norm(a1)
            v[6:9] = a3 / np.linalg.norm(a3)
        out.append(' '.join('%.15g' % x for x in v) + '\n')
    with open(a.out, 'w') as f:
        f.writelines(out)


if __name__ == '__main__':
    main()
