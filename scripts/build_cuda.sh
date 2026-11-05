#!/usr/bin/env bash
# Build CUDA versions of both codes (for a cluster node with an NVIDIA GPU and nvcc).
#
#   scripts/build_cuda.sh [COMMIT]
#
# The original anm-oxdna CUDA code uses texture references, which were removed in CUDA 12: build it
# with a CUDA 11.x toolkit (or skip it and compare the new CUDA build against the old CPU build,
# which is the default in cuda_tests.md).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
COMMIT="${1:-${NEW_COMMIT:-f802bea4}}"
"$REPO/scripts/build_new.sh" "$COMMIT" "$REPO/build/new_cuda_$COMMIT" -DCUDA=ON
if [ "${BUILD_OLD_CUDA:-0}" = "1" ]; then
    "$REPO/scripts/build_old.sh" "${OLD_SRC:-$HOME/software/anm-oxdna/oxDNA}" "$REPO/build/old_cuda" -DCUDA=ON -DCUDA_COMMON_ARCH=OFF
fi
