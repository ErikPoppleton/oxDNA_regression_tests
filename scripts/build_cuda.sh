#!/usr/bin/env bash
# Build the current oxDNA for the CUDA tests: CPU (double), CUDA mixed precision and CUDA double precision,
# out of tree, from a checkout that is left untouched.
#
#   scripts/build_cuda.sh <oxDNA checkout> <build root>       -> <build root>/{cpu,cuda_mixed,cuda_double}/bin/oxDNA
#
# scripts/cluster_cuda.sh does the same and then runs the tests; use this one to build by hand.
# Extra CMake arguments can be given with CMAKE_EXTRA (default -DCUDA_COMMON_ARCH=OFF: build for the GPU
# of the machine).
set -euo pipefail
SRC="$(cd "$1" && pwd)"
ROOT="$(mkdir -p "$2" && cd "$2" && pwd)"
JOBS="${JOBS:-8}"
CMAKE_EXTRA="${CMAKE_EXTRA:--DCUDA_COMMON_ARCH=OFF}"
# shellcheck disable=SC2086
cmake -S "$SRC" -B "$ROOT/cpu" -DCUDA=OFF $CMAKE_EXTRA && make -C "$ROOT/cpu" -j"$JOBS" oxDNA
# shellcheck disable=SC2086
cmake -S "$SRC" -B "$ROOT/cuda_mixed" -DCUDA=ON $CMAKE_EXTRA && make -C "$ROOT/cuda_mixed" -j"$JOBS" oxDNA
# shellcheck disable=SC2086
cmake -S "$SRC" -B "$ROOT/cuda_double" -DCUDA=ON -DCUDA_DOUBLE=ON $CMAKE_EXTRA && make -C "$ROOT/cuda_double" -j"$JOBS" oxDNA
echo "built: $ROOT/cpu/bin/oxDNA $ROOT/cuda_mixed/bin/oxDNA $ROOT/cuda_double/bin/oxDNA"
