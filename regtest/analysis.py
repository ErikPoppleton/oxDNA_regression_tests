"""Statistics used by the regression tests: block averages, alignment and RMSF."""
import numpy as np

SIM_LENGTH_A = 8.518  # Angstrom per simulation length unit (ANM convention, ANMUtils/models.py)


def block_average(x, n_blocks=10):
    """Mean and standard error of a correlated time series from n_blocks non-overlapping blocks."""
    x = np.asarray(x, dtype=float)
    n = len(x) // n_blocks
    if n < 1:
        return float(np.mean(x)), float('nan')
    blocks = x[:n * n_blocks].reshape(n_blocks, n).mean(axis=1)
    return float(np.mean(x[:n * n_blocks])), float(blocks.std(ddof=1) / np.sqrt(n_blocks))


def unwrap_to_reference(frames, box, ref_index=0):
    """Make every particle the minimum image of itself in the previous frame (removes PBC jumps and
    the per-strand shifts oxDNA applies when fix_diffusion is on), then the minimum image with respect
    to particle ref_index inside each frame is NOT applied: molecules spanning strands are kept whole by
    chaining each frame to the previous one, starting from a whole first frame."""
    out = frames.copy()
    for t in range(1, len(out)):
        d = out[t] - out[t - 1]
        d -= box[t] * np.rint(d / box[t])
        out[t] = out[t - 1] + d
    return out


def make_whole(coords, box):
    """Chain particles i -> i-1 by minimum image so a molecule split by periodic boundaries is whole."""
    c = coords.copy()
    for i in range(1, len(c)):
        d = c[i] - c[i - 1]
        d -= box * np.rint(d / box)
        c[i] = c[i - 1] + d
    return c


def kabsch(P, Q):
    """Rotation R minimising |(P - <P>) R - (Q - <Q>)|."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    return U @ D @ Vt


def align(frames, ref, n_iter=3):
    """Align all frames on ref (iterating the reference to the mean structure); returns aligned frames
    and the final mean structure."""
    ref = ref - ref.mean(0)
    for _ in range(n_iter):
        aligned = np.empty_like(frames)
        for t, X in enumerate(frames):
            R = kabsch(X, ref)
            aligned[t] = (X - X.mean(0)) @ R
        ref = aligned.mean(0)
    return aligned, ref


def rmsf(frames, n_blocks=5):
    """Per-particle RMSF (input units) of pre-aligned frames, plus standard error from n_blocks blocks."""
    mean = frames.mean(0)
    full = np.sqrt(((frames - mean) ** 2).sum(-1).mean(0))
    n = len(frames) // n_blocks
    if n < 2:
        return full, np.full_like(full, np.nan)
    per_block = []
    for b in range(n_blocks):
        fb = frames[b * n:(b + 1) * n]
        per_block.append(np.sqrt(((fb - mean) ** 2).sum(-1).mean(0)))
    per_block = np.array(per_block)
    return full, per_block.std(0, ddof=1) / np.sqrt(n_blocks)


def trajectory_rmsf(traj_data, boxes, idx, equil_frac=0.2, n_blocks=5):
    """RMSF (Angstrom) of particles idx, aligned on idx. traj_data: (n_frames, N, >=3)."""
    pos = traj_data[:, :, :3]
    start = int(len(pos) * equil_frac)
    pos = unwrap_to_reference(pos, boxes)
    pos[0] = make_whole(pos[0], boxes[0])
    pos = unwrap_to_reference(pos, boxes)  # re-chain from the now-whole first frame
    sel = pos[start:, idx, :]
    aligned, _ = align(sel, sel[0])
    r, se = rmsf(aligned, n_blocks=n_blocks)
    return r * SIM_LENGTH_A, se * SIM_LENGTH_A
