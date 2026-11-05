# Running the ANM regression tests on CUDA

These steps need a node with an NVIDIA GPU, `nvcc`, CMake >= 3.18 and Python 3 with numpy.
Nothing here has been run yet: the CPU suite was developed on a laptop without a GPU.

## 1. Build

```bash
git clone https://github.com/ErikPoppleton/oxDNA_regression_tests.git
cd oxDNA_regression_tests

# reference: the original code, CPU, double precision (the default reference for CUDA runs)
scripts/build_old.sh /path/to/anm-oxdna/oxDNA build/old

# code under test, CUDA (mixed precision). NEW_REPO is the current oxDNA checkout;
# the build happens in a detached worktree, so the checkout is not touched.
NEW_REPO=/path/to/oxDNA scripts/build_new.sh f802bea4 build/new_cuda -DCUDA=ON

# optional: the original code on CUDA. It uses CUDA texture references, which were removed in
# CUDA 12, so this needs a CUDA 11.x toolkit.
scripts/build_old.sh /path/to/anm-oxdna/oxDNA build/old_cuda -DCUDA=ON -DCUDA_COMMON_ARCH=OFF
```

Replace `f802bea4` with the commit you want to test (e.g. the head of `anm_refactor`).

## 2. Run

New CUDA against old CPU (recommended; the reference is the double-precision original):

```bash
./run_tests.py --old build/old/bin/oxDNA \
               --new build/new_cuda/build/bin/oxDNA --backend CUDA \
               --jobs 2 --out results/CUDA
```

Quick subset (a few minutes):

```bash
./run_tests.py --old build/old/bin/oxDNA --new build/new_cuda/build/bin/oxDNA --backend CUDA --quick --jobs 2
```

New CUDA against old CUDA (only if the old CUDA build worked):

```bash
./run_tests.py --old build/old_cuda/bin/oxDNA --old-backend CUDA \
               --new build/new_cuda/build/bin/oxDNA --backend CUDA --jobs 2
```

New CUDA against new CPU (the same code on both backends; useful to separate GPU problems from
model differences):

```bash
./run_tests.py --old build/new_f802bea4/build/bin/oxDNA \
               --new build/new_cuda/build/bin/oxDNA --backend CUDA --jobs 2
```

(in this mode the `old` inputs of each case are still used for the reference runs, so cases whose
old/new inputs differ, i.e. `*/ANMT` which runs `ACT` on the old side, will fail with an unknown
interaction; restrict with `--only '*/DNANM' '*/AC'` etc.)

`--jobs` is the number of simultaneous simulations; with one GPU keep it small (2-4) so the
GPU runs do not time out. Results go to `results/CUDA/summary.md` (table) and `results.json`.

Please send back `results/CUDA/summary.md` and, for failures, the `stdout.txt`/`log.dat` of the
failing run directory (`results/CUDA/<system>/<case>/<stage>/<old|new>/`).

## 3. What is compared on CUDA

- `energy0`: step-0 energy of the perturbed configuration. Tolerance relaxed to 1e-4 relative
  (mixed precision).
- `det`: skipped (trajectories from different backends/RNGs cannot be compared step by step).
- `stat`: block-averaged potential energy and per-residue RMSF, as on the CPU.

## 4. Expected status at f802bea4

| case | expected on CUDA | why |
|---|---|---|
| `*/AC` (ubiquitin, cln025, kdpg) | XFAIL (error) | there is no CUDA AC interaction; protein-only ANM on the GPU goes through DNANM |
| `*/ANMT`, `*/DNACT` (all systems) | XFAIL (error) | `CUDADNANMInteraction` throws for DNACT: bending/torsion are not implemented on the GPU |
| `*/RNANM`, `*/RNACT` (ubiquitin, cln025, ubq_rna) | XFAIL (error) | no CUDA RNANM/RNACT in the current code |
| `ubiquitin/DNANM`, `kdpg/DNANM` | PASS | protein-only DNANM is supported on CUDA |
| `ubq_dna/DNANM`, `cage/DNANM` | PASS | protein + DNA on CUDA |

Anything marked PASS above that fails is a real CUDA problem. Cases marked XFAIL that start to
pass will be reported as XPASS (not an error): update the `xfail` entries in the system's
`config.json` when CUDA support is added.

Notes on the original CUDA code, if you use it as the reference: its `CUDAInteractionFactory` maps
`RNANM` to the **DNA** GPU interaction (`CUDADNANMInteraction`), so old-CUDA RNANM results are not a
valid reference for RNANM; and it has no protein-only `AC`/`ACT` GPU interaction.
