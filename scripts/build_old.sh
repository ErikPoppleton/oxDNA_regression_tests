#!/usr/bin/env bash
# Build the original ANM implementation (anm-oxdna fork of an old oxDNA) out of tree.
#
#   scripts/build_old.sh [SRC] [BUILD_DIR] [extra cmake args...]
#
#   SRC        anm-oxdna/oxDNA checkout       (default: $OLD_SRC or ~/software/anm-oxdna/oxDNA)
#   BUILD_DIR  where to build                 (default: <repo>/build/old)
#   e.g. CUDA: scripts/build_old.sh ~/anm-oxdna/oxDNA build/old_cuda -DCUDA=ON -DCUDA_COMMON_ARCH=OFF
#
# The source tree is never modified. Precision is chosen at run time in the old code
# (backend_precision = double for CPU, mixed for CUDA), so no precision flag is needed here.
#
# Build fixes: none were needed on macOS arm64 (Apple clang 21, CMake 3.31) except
# -DCMAKE_POLICY_VERSION_MINIMUM=3.5, required by CMake >= 4 because the old CMakeLists.txt asks
# for cmake_minimum_required(2.6). If a patch is ever needed it goes into patches/ and is applied to
# a copy of the sources (see APPLY_PATCHES below).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${1:-${OLD_SRC:-$HOME/software/anm-oxdna/oxDNA}}"
BUILD="${2:-$REPO/build/old}"
shift $(( $# > 2 ? 2 : $# )) || true
JOBS="${JOBS:-8}"

SRC="$(cd "$SRC" && pwd)"
mkdir -p "$BUILD"
BUILD="$(cd "$BUILD" && pwd)"

# optional patches (none needed at the moment): copy the sources and apply every patches/anm-oxdna-*.patch
if [ "${APPLY_PATCHES:-0}" = "1" ] && ls "$REPO"/patches/anm-oxdna-*.patch >/dev/null 2>&1; then
    rm -rf "$BUILD/src_copy"
    cp -R "$SRC" "$BUILD/src_copy"
    for p in "$REPO"/patches/anm-oxdna-*.patch; do
        echo "applying $p"
        patch -d "$BUILD/src_copy" -p1 < "$p"
    done
    SRC="$BUILD/src_copy"
fi

cmake -S "$SRC" -B "$BUILD" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 "$@"
make -C "$BUILD" -j"$JOBS" oxDNA DNAnalysis
echo "old binary: $BUILD/bin/oxDNA"
