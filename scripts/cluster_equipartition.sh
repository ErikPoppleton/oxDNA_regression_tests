#!/usr/bin/env bash
#SBATCH --job-name=oxdna_equipartition
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=oxdna_equipartition_%j.log
#
# Equipartition check of the CUDA MD backend for protein systems (tools/equipartition.py): measures the
# translational and rotational temperatures of amino acids and nucleotides from the velocities in the
# trajectories, on the CPU and CUDA (mixed, double) builds of the branch and on the CUDA DNANM of master.
#
#   sbatch scripts/cluster_equipartition.sh <oxDNA checkout> <outdir>     # submit from oxDNA_regression_tests
#
# <oxDNA checkout>  the oxDNA source tree on the branch to test (anm_refactor); built out of tree, not modified.
# <outdir>          builds, runs, equipartition.md and a tarball go here.
#
# Environment variables (all optional):
#   MASTER_SRC       oxDNA source tree of master. Default: exported from the "master" ref of <oxDNA checkout>
#                    (git archive, the checkout is not touched). MASTER_SRC=none skips the master runs.
#   CLUSTER_MODULES  modules to load before building (as in cluster_cuda.sh)
#   PYTHON           python3 with numpy (default: python3)
#   JOBS             parallel build jobs / CPU simulations (default: $SLURM_CPUS_PER_TASK or nproc)
#   GPU_JOBS         simultaneous simulations on the GPU (default: 4)
#   CMAKE_EXTRA      extra CMake arguments (default -DCUDA_COMMON_ARCH=OFF: build for the GPU of the node)
#   SKIP_BUILD=1     reuse the builds already in <outdir>/build
#   ONLY             space separated run names of tools/equipartition.py (default: all)
set -uo pipefail

if [ $# -lt 2 ]; then
    sed -n '2,26p' "$0"
    exit 2
fi
SRC="$(cd "$1" && pwd)"
OUT="$(mkdir -p "$2" && cd "$2" && pwd)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -f "$REPO/tools/equipartition.py" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
    [ -f "$REPO/tools/equipartition.py" ] || REPO="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
[ -f "$REPO/tools/equipartition.py" ] || { echo "cannot find tools/equipartition.py (submit from the oxDNA_regression_tests directory)"; exit 2; }
PYTHON="${PYTHON:-python3}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 8)}}"
GPU_JOBS="${GPU_JOBS:-4}"
CMAKE_EXTRA="${CMAKE_EXTRA:--DCUDA_COMMON_ARCH=OFF}"
MASTER_SRC="${MASTER_SRC:-}"
LOG="$OUT/cluster_equipartition.log"
SUMMARY="$OUT/SUMMARY.md"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

if [ -n "${CLUSTER_MODULES:-}" ]; then
    type module > /dev/null 2>&1 || source /etc/profile.d/modules.sh 2> /dev/null || source /usr/share/lmod/lmod/init/bash 2> /dev/null
    # shellcheck disable=SC2086
    module load $CLUSTER_MODULES || { echo "module load $CLUSTER_MODULES failed"; exit 2; }
fi
"$PYTHON" -c "import numpy" || { echo "PYTHON=$PYTHON has no numpy"; exit 2; }

: > "$LOG"
if [ "$MASTER_SRC" != "none" ] && [ -z "$MASTER_SRC" ]; then
    MASTER_SRC="$OUT/master_src"
    if [ ! -f "$MASTER_SRC/CMakeLists.txt" ]; then
        log "exporting master from $SRC"
        mkdir -p "$MASTER_SRC"
        git -C "$SRC" archive master | tar -x -C "$MASTER_SRC" || { log "git archive master failed: skipping master"; MASTER_SRC=none; }
    fi
fi
MASTER_REV="$( [ "$MASTER_SRC" != none ] && git -C "$SRC" rev-parse --short master 2>/dev/null)"

{
    echo "# oxDNA equipartition check"
    echo
    echo "- date: $(date)"
    echo "- host: $(hostname)"
    echo "- branch source: $SRC ($(git -C "$SRC" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "$SRC" rev-parse --short HEAD 2>/dev/null))"
    echo "- master source: $MASTER_SRC ${MASTER_REV:+(master $MASTER_REV)}"
    echo "- test repo: $REPO ($(git -C "$REPO" rev-parse --short HEAD 2>/dev/null))"
    echo "- GPU: $(nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader 2>/dev/null | head -1)"
    echo "- nvcc: $(nvcc --version 2>/dev/null | tail -1)"
    echo
} > "$SUMMARY"

build() {  # name, source, cmake args...
    local name="$1" src="$2"; shift 2
    local dir="$OUT/build/$name"
    if [ "${SKIP_BUILD:-0}" = "1" ] && [ -x "$dir/bin/oxDNA" ]; then
        log "reusing build $name"
        return 0
    fi
    log "building $name ($*)"
    mkdir -p "$dir"
    # shellcheck disable=SC2086
    if cmake -S "$src" -B "$dir" "$@" $CMAKE_EXTRA > "$dir/cmake.log" 2>&1 && make -C "$dir" -j"$JOBS" oxDNA > "$dir/make.log" 2>&1; then
        log "built $dir/bin/oxDNA"
    else
        log "BUILD FAILED: $name (see $dir/make.log)"
        echo "- **build $name FAILED**" >> "$SUMMARY"
        tail -30 "$dir/make.log" >> "$LOG" 2>/dev/null
        return 1
    fi
}

build cpu "$SRC" -DCUDA=OFF || exit 1
build cuda_mixed "$SRC" -DCUDA=ON || exit 1
build cuda_double "$SRC" -DCUDA=ON -DCUDA_DOUBLE=ON
ARGS=(--out "$OUT/runs" --jobs "$JOBS" --gpu-jobs "$GPU_JOBS"
      --branch-cpu "$OUT/build/cpu/bin/oxDNA" --branch-mixed "$OUT/build/cuda_mixed/bin/oxDNA")
[ -x "$OUT/build/cuda_double/bin/oxDNA" ] && ARGS+=(--branch-double "$OUT/build/cuda_double/bin/oxDNA")
if [ "$MASTER_SRC" != "none" ]; then
    build master_mixed "$MASTER_SRC" -DCUDA=ON && ARGS+=(--master-mixed "$OUT/build/master_mixed/bin/oxDNA")
fi
# shellcheck disable=SC2206
[ -n "${ONLY:-}" ] && ARGS+=(--only $ONLY)

log "running: equipartition.py run ${ARGS[*]}"
"$PYTHON" "$REPO/tools/equipartition.py" run "${ARGS[@]}" 2>&1 | tee -a "$LOG"
RUN_RC=${PIPESTATUS[0]}
log "runs finished (exit code $RUN_RC); analysing"
"$PYTHON" "$REPO/tools/equipartition.py" analyze --out "$OUT/runs" > /dev/null 2>> "$LOG"
cat "$OUT/runs/equipartition.md" >> "$SUMMARY" 2>/dev/null
[ "$RUN_RC" -ne 0 ] && echo "**$(grep -c '^FAILED' "$LOG") run(s) failed, see cluster_equipartition.log**" >> "$SUMMARY"

# tarball: everything but the trajectories (per-particle averages are in per_particle.npz) and the builds
TAR="$OUT/oxdna_equipartition_$(date +%Y%m%d_%H%M).tar.gz"
(cd "$OUT" && find . \( -path ./build -o -path ./master_src \) -prune -o -type f \( -name '*.md' -o -name '*.json' \
    -o -name '*.npz' -o -name 'input' -o -name 'energy.dat' -o -name 'stdout.txt' -o -name '*.log' \) -print \
    | tar -czf "$TAR" -T - 2>/dev/null; tar -czf build_logs.tar.gz build/*/cmake.log build/*/make.log 2>/dev/null)
log "summary: $SUMMARY"
log "tarball: $TAR"
cat "$SUMMARY"
exit "$RUN_RC"
