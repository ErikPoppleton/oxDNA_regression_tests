"""Minimal readers for oxDNA files (no dependency on either oxDNA code)."""
import numpy as np


def read_energy(path, with_gpu=False):
    """Return the energy file as an array (columns: time, potential, kinetic, total, ...).
    With with_gpu, also return the 'CUDA_energy:' values that CUDA_print_energy appends (nan if absent)."""
    rows, gpu = [], []
    with open(path) as f:
        for line in f:
            g = float('nan')
            if 'CUDA_energy:' in line:
                line, rest = line.split('CUDA_energy:', 1)
                try:
                    g = float(rest.split()[0])
                except (ValueError, IndexError):
                    pass
            s = line.split()
            if not s or s[0].startswith('#'):
                continue
            try:
                rows.append([float(x) for x in s])
                gpu.append(g)
            except ValueError:
                continue
    if not rows:
        arr = np.zeros((0, 4))
    else:
        n = min(len(r) for r in rows)
        arr = np.array([r[:n] for r in rows])
    return (arr, np.array(gpu)) if with_gpu else arr


def read_columns(path):
    """Plain whitespace-separated numeric file (e.g. a data_output observable file)."""
    rows = []
    with open(path) as f:
        for line in f:
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    return np.array(rows)


def read_trajectory(path, cols=9):
    """Return (times, boxes, data) with data of shape (n_frames, N, cols)."""
    times, boxes, frames, cur = [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('t ='):
                if cur:
                    frames.append(cur)
                cur = []
                times.append(float(line.split()[2]))
            elif line.startswith('b ='):
                boxes.append([float(x) for x in line.split()[2:5]])
            elif line.startswith('E ='):
                continue
            elif line.strip():
                cur.append([float(x) for x in line.split()[:cols]])
    if cur:
        frames.append(cur)
    n = min(len(fr) for fr in frames) if frames else 0
    frames = [fr for fr in frames if len(fr) == n]  # drop a truncated last frame, if any
    return np.array(times[:len(frames)]), np.array(boxes[:len(frames)]), np.array(frames)


def read_conf(path):
    t, b, d = read_trajectory(path, cols=15)
    return b[0], d[0]


def read_par(path):
    """Return the springs [(i, j)] of a .par file."""
    springs = []
    with open(path) as f:
        f.readline()
        for line in f:
            s = line.split()
            if len(s) >= 5:
                springs.append((int(s[0]), int(s[1])))
    return springs


def read_old_topology(path):
    """Return (header, particles) for an old-style topology; particles = (strand, code, n3, n5, extra ints)."""
    with open(path) as f:
        header = f.readline().split()
        parts = []
        for line in f:
            s = line.split()
            if not s or s[0].startswith('#'):
                continue
            parts.append((int(s[0]), s[1], int(s[2]), int(s[3]), [int(x) for x in s[4:]]))
    return header, parts
