# CUDA tests (for a GPU node)

This tests the CUDA implementation of the protein models (`DNANM`, `DNACT`, `RNANM`, `RNACT`, `NANM`,
`NACT`; `src/CUDA/Interactions/CUDAANMInteraction.cu`) and of oxNA (`NA`) against the CPU code of the
**same** oxDNA build. Both CUDA precisions are tested. None of this has run on a GPU yet: it was written
without one.

## What to run

On the cluster (the test repo has to be there: push it from the laptop and clone it, or rsync it):

```bash
git clone https://github.com/ErikPoppleton/oxDNA_regression_tests.git
cd oxDNA_regression_tests
git -C /path/to/oxDNA fetch && git -C /path/to/oxDNA checkout anm_refactor

# module names are site specific: list what your cluster needs (CUDA, CMake >= 3.18, a C++14 compiler)
export CLUSTER_MODULES="cuda/12.4 cmake/3.29 gcc/12.3"
sbatch -p <gpu partition> scripts/cluster_cuda.sh /path/to/oxDNA $SCRATCH/oxdna_anm_cuda
```

Submit from the `oxDNA_regression_tests` directory (the script finds `run_tests.py` through
`$SLURM_SUBMIT_DIR`). On an interactive GPU node, `bash scripts/cluster_cuda.sh /path/to/oxDNA <outdir>` does the same.
Other options (`PYTHON`, `JOBS`, `GPU_JOBS`, `CMAKE_EXTRA`, `SKIP_BUILD`, `SKIP_DOUBLE`,
`CONTINUE_ON_SMOKE_FAILURE`, `OLD_OXDNA`) are documented at the top of the script. Python 3 with numpy is needed.

The script:

1. builds `<outdir>/build/cpu` (CPU, double) and `<outdir>/build/cuda_mixed` (`-DCUDA=ON`);
2. runs the **smoke test** (`--smoke`, 8 cases: step-0 energies and a 2000-step NVE run, 1-2 min). If it
   fails the script stops here, so a broken build is reported in minutes;
3. builds `<outdir>/build/cuda_double` (`-DCUDA=ON -DCUDA_DOUBLE=ON`);
4. runs the full suite in mixed precision and in double precision (CPU build as the reference, `--same-code`);
5. runs the oxDNA repository's own `test/CUDA/NA/test_cuda_na.py` with both builds;
6. writes `<outdir>/SUMMARY.md` and `<outdir>/oxdna_anm_cuda_results_<date>.tar.gz`.

Please send back the tarball (it includes `SUMMARY.md`, every `summary.md`/`results.json`, and the input, log
and energy files of every run) and `<outdir>/build_logs.tar.gz` if a build failed.

Expected time with 8 CPU cores and one GPU: builds 20-40 min, smoke 1-2 min, each full suite about
15-25 min (the CPU reference runs dominate).

Individual commands, once the binaries exist:

```bash
B=$SCRATCH/oxdna_anm_cuda/build
./run_tests.py --same-code --old $B/cpu/bin/oxDNA --new $B/cuda_mixed/bin/oxDNA  --backend CUDA --cuda-precision mixed  --smoke
./run_tests.py --same-code --old $B/cpu/bin/oxDNA --new $B/cuda_mixed/bin/oxDNA  --backend CUDA --cuda-precision mixed
./run_tests.py --same-code --old $B/cpu/bin/oxDNA --new $B/cuda_double/bin/oxDNA --backend CUDA --cuda-precision double
./run_tests.py --same-code ... --only 'ubq_hybrid/*' 'na_*'      # a subset
./run_tests.py --old anm-oxdna/bin/oxDNA --new $B/cuda_mixed/bin/oxDNA --backend CUDA   # original CPU code vs current CUDA
```

## What is compared on CUDA

In `--same-code` mode both sides use the same input (the `new` inputs of a case, plus its `new@CUDA`
substitutions). The only difference is `backend = CPU` and `backend = CUDA`, with `use_edge = true`,
`CUDA_sort_every = 0` and `CUDA_list = verlet`; double precision uses `use_edge = false`.
Protein-only `AC`/`ACT` cases run as `DNANM`/`DNACT` with the protein-only topology, which the CPU
suite shows is the same simulation. Duplicate cases (ubiquitin `AC`/`ANMT`/`ACT`, kdpg `AC`/`ACT`) are
not run in this mode.

| check | what it means on CUDA | passes if |
|---|---|---|
| `energy0` | perturbed configuration, external forces off. One step with `CUDA_print_energy = true`: the energy line at t = dt contains the host (CPU-code) energy and the energy summed on the GPU, for the same configuration. This is a **direct test of the GPU energies**. The host energy at t = 0 is also compared with the CPU run, which is trivially equal. | GPU - host < 2e-5 + 1e-4 \|E\| per particle (the GPU uses analytic functions where the CPU uses interpolation meshes: ~1e-5) |
| `nve` | 20 000 steps, dt = 0.001, `thermostat = no`, from `equilibrated.dat`. The energy file is computed by the CPU code from the GPU positions, so a drift of U + K (+ half the `force_energy` of the mutual-trap pairs) means **the GPU forces/torques are not the gradient of the energy**. Drift = linear fit over the run; fluctuations = rms about the fit. | drift_GPU <= 3 x drift_CPU + 2e-4 and fluct_GPU <= 3 x fluct_CPU + 2e-4 (per particle) |
| `stat` | 4-8 independent replicas per backend: mean potential energy and per-residue protein RMSF (aligned) | as on the CPU (README): \|dE\| <= 4 sigma (4.5 for nucleic systems), RMSF per residue / mean within 4 sigma |
| `det` | skipped: trajectories from different backends cannot be compared step by step | - |

Typical CPU values (laptop): NVE drift 1e-5 to 4e-4 per particle over 20 time units, fluctuations 1e-5 to 6e-5.

## Expected status at anm_refactor ca0afa94

| precision | cases | expected |
|---|---|---|
| mixed | all ANM cases (DNANM, DNACT, RNANM, RNACT, NANM, NACT; protein only, DNA/RNA hybrids, cage, DNA+RNA+protein) and all oxNA (`na_*`) cases | PASS |
| double | `na_*` (oxNA) | PASS |

Statuses: PASS / FAIL (the comparison failed) / ERROR (a run crashed: see its `log.dat`) / XFAIL (expected failure) /
XPASS (expected to fail but passed: e.g. once ANM works in double precision, remove the `*@CUDA_double`
entries from the `xfail` blocks of the `config.json` files) / SKIP.

## Reading the output

- `SUMMARY.md`: one section per suite with the totals, plus the rows that are FAIL, ERROR or XPASS.
- `results/<suite>/summary.md`: every check of every case. `ref` is the CPU build and `test` is the CUDA build.
- `results/<suite>/results.json`: all numbers, including per-replica energies, per-residue RMSFs and NVE drifts.
- `results/<suite>/<system>/<case>/<stage>/{old,new}[/repN]/`: run directories (`input`, `log.dat`,
  `energy.dat`, `ext.dat`). `old` is the CPU reference and `new` is CUDA.
- `test_cuda_na_{mixed,double}.log`: output of oxDNA's own NA CUDA tests.
