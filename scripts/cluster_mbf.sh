#!/usr/bin/env bash
#SBATCH --job-name=oxdna_mbf_cuda
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=oxdna_mbf_cuda_%j.log
#
# CPU/CUDA consistency of the capped backbone potential (max_backbone_force, max_backbone_rate): see
# tools/mbf_cuda_check.py. Builds the CPU, CUDA mixed and CUDA double versions of oxDNA and the ANNaMo plugins (CPU and
# CUDA), runs the checks and packs the results.
#
#   sbatch scripts/cluster_mbf.sh <oxDNA checkout> <outdir>     # submit from oxDNA_regression_tests
#
# <oxDNA checkout>  the oxDNA source tree on the branch to test (max_backbone_rate). It is built out of tree, but the
#                   ANNaMo plugins are written to its contrib/tostiguerra folder (as for any oxDNA plugin build) and
#                   then copied to <outdir>/plugins_*.
# <outdir>          builds, runs, mbf_cuda_check.md and a tarball go here.
#
# Environment variables (all optional):
#   CLUSTER_MODULES  modules to load before building (as in cluster_cuda.sh)
#   PYTHON           python3 with numpy (default: python3)
#   JOBS             parallel build jobs (default: $SLURM_CPUS_PER_TASK or nproc)
#   CMAKE_EXTRA      extra CMake arguments (default -DCUDA_COMMON_ARCH=OFF: build for the GPU of the node)
#   SKIP_BUILD=1     reuse the builds already in <outdir>/build
#   SKIP_DOUBLE=1    do not build/test the CUDA_DOUBLE version
#   ONLY             space separated case names of tools/mbf_cuda_check.py (default: all)
set -uo pipefail

if [ $# -lt 2 ]; then
    sed -n '2,25p' "$0"
    exit 2
fi
SRC="$(cd "$1" && pwd)"
OUT="$(mkdir -p "$2" && cd "$2" && pwd)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -f "$REPO/tools/mbf_cuda_check.py" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
    [ -f "$REPO/tools/mbf_cuda_check.py" ] || REPO="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
[ -f "$REPO/tools/mbf_cuda_check.py" ] || { echo "cannot find tools/mbf_cuda_check.py (submit from the oxDNA_regression_tests directory)"; exit 2; }
PYTHON="${PYTHON:-python3}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 8)}}"
CMAKE_EXTRA="${CMAKE_EXTRA:--DCUDA_COMMON_ARCH=OFF}"
LOG="$OUT/cluster_mbf.log"
SUMMARY="$OUT/SUMMARY.md"
PLUGIN_DIR="$SRC/contrib/tostiguerra"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

if [ -n "${CLUSTER_MODULES:-}" ]; then
    type module > /dev/null 2>&1 || source /etc/profile.d/modules.sh 2> /dev/null || source /usr/share/lmod/lmod/init/bash 2> /dev/null
    # shellcheck disable=SC2086
    module load $CLUSTER_MODULES || { echo "module load $CLUSTER_MODULES failed"; exit 2; }
fi
"$PYTHON" -c "import numpy" || { echo "PYTHON=$PYTHON has no numpy"; exit 2; }

: > "$LOG"
{
    echo "# oxDNA max_backbone_force CPU/CUDA check"
    echo
    echo "- date: $(date)"
    echo "- host: $(hostname)"
    echo "- oxDNA source: $SRC ($(git -C "$SRC" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "$SRC" rev-parse --short HEAD 2>/dev/null))"
    echo "- test repo: $REPO ($(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null))"
    echo "- GPU: $(nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader 2>/dev/null | head -1)"
    echo "- nvcc: $(nvcc --version 2>/dev/null | tail -1)"
    echo
} > "$SUMMARY"

build() {  # name, cmake args...
    local name="$1"; shift
    local dir="$OUT/build/$name"
    if [ "${SKIP_BUILD:-0}" = "1" ] && [ -x "$dir/bin/oxDNA" ]; then
        log "reusing build $name"
        return 0
    fi
    log "building $name ($*)"
    mkdir -p "$dir"
    # shellcheck disable=SC2086
    if cmake -S "$SRC" -B "$dir" "$@" $CMAKE_EXTRA > "$dir/cmake.log" 2>&1 && make -C "$dir" -j"$JOBS" oxDNA > "$dir/make.log" 2>&1; then
        log "built $dir/bin/oxDNA"
    else
        log "BUILD FAILED: $name (see $dir/make.log)"
        echo "- **build $name FAILED**" >> "$SUMMARY"
        tail -30 "$dir/make.log" >> "$LOG" 2>/dev/null
        return 1
    fi
}

plugins() {  # build name, destination, targets...
    local name="$1" dest="$2"; shift 2
    mkdir -p "$dest"
    for target in "$@"; do
        if make -C "$OUT/build/$name" "$target" >> "$OUT/build/$name/plugins.log" 2>&1 && [ -f "$PLUGIN_DIR/$target.so" ]; then
            cp "$PLUGIN_DIR/$target.so" "$dest/"
            log "built plugin $target ($name)"
        else
            log "PLUGIN BUILD FAILED: $target ($name), see build/$name/plugins.log"
            echo "- **plugin $target ($name) FAILED** (see build/$name/plugins.log)" >> "$SUMMARY"
        fi
    done
}

build cpu -DCUDA=OFF || exit 1
plugins cpu "$OUT/plugins_cpu" ANNaMoInteraction
build cuda_mixed -DCUDA=ON || exit 1
plugins cuda_mixed "$OUT/plugins_cuda" ANNaMoInteraction CUDAANNaMoInteraction
ARGS=(--out "$OUT/runs" --cpu "$OUT/build/cpu/bin/oxDNA" --mixed "$OUT/build/cuda_mixed/bin/oxDNA"
      --cpu-plugins "$OUT/plugins_cpu" --cuda-plugins "$OUT/plugins_cuda")
if [ "${SKIP_DOUBLE:-0}" != "1" ]; then
    build cuda_double -DCUDA=ON -DCUDA_DOUBLE=ON && ARGS+=(--double "$OUT/build/cuda_double/bin/oxDNA")
fi
# shellcheck disable=SC2206
[ -n "${ONLY:-}" ] && ARGS+=(--only $ONLY)

log "running: mbf_cuda_check.py run ${ARGS[*]}"
"$PYTHON" "$REPO/tools/mbf_cuda_check.py" run "${ARGS[@]}" 2>&1 | tee -a "$LOG"
"$PYTHON" "$REPO/tools/mbf_cuda_check.py" analyze --out "$OUT/runs" > /dev/null 2>> "$LOG"
RC=$?
cat "$OUT/runs/mbf_cuda_check.md" >> "$SUMMARY" 2>/dev/null
grep -q '^FAILED' "$LOG" && echo "**$(grep -c '^FAILED' "$LOG") run(s) failed, see cluster_mbf.log**" >> "$SUMMARY"

# tarball: everything but the trajectories and the builds
TAR="$OUT/oxdna_mbf_cuda_$(date +%Y%m%d_%H%M).tar.gz"
(cd "$OUT" && find . \( -path ./build -o -path ./plugins_cpu -o -path ./plugins_cuda \) -prune -o -type f \( -name '*.md' \
    -o -name '*.json' -o -name 'input' -o -name 'energy.dat' -o -name 'log.dat' -o -name 'stdout.txt' -o -name '*.log' \) -print \
    | tar -czf "$TAR" -T - 2>/dev/null; tar -czf build_logs.tar.gz build/*/cmake.log build/*/make.log build/*/plugins.log 2>/dev/null)
log "summary: $SUMMARY"
log "tarball: $TAR"
cat "$SUMMARY"
exit "$RC"
