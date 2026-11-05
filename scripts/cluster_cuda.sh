#!/usr/bin/env bash
#SBATCH --job-name=oxdna_anm_cuda_tests
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=oxdna_anm_cuda_tests_%j.log
#
# Build oxDNA (CPU, CUDA mixed, CUDA double) and run the CUDA regression suite on a GPU node.
#
#   sbatch scripts/cluster_cuda.sh <oxDNA checkout> <outdir>        # or: bash scripts/cluster_cuda.sh ... on a GPU node
#
# <oxDNA checkout>  the oxDNA source tree to test (check out anm_refactor first). It is built out of tree
#                   in <outdir>/build/*; the checkout itself is not modified.
# <outdir>          everything goes here: builds, results, SUMMARY.md and a tarball.
#
# Environment variables (all optional):
#   CLUSTER_MODULES  modules to load before building, e.g. "cuda/12.4 cmake/3.29 gcc/12.3". Module names
#                    are site specific, so none are hard-coded; if unset, nothing is loaded and nvcc/cmake/
#                    a C++ compiler must already be in PATH.
#   PYTHON           python3 with numpy (default: python3)
#   JOBS             parallel build jobs / parallel CPU simulations (default: $SLURM_CPUS_PER_TASK or nproc)
#   GPU_JOBS         simultaneous simulations on the GPU (default: 4)
#   CMAKE_EXTRA      extra CMake arguments for all builds (e.g. -DCUDA_COMMON_ARCH=ON to build for all archs;
#                    the default -DCUDA_COMMON_ARCH=OFF compiles for the GPU of the node)
#   SKIP_BUILD=1     reuse the builds already in <outdir>/build
#   SKIP_DOUBLE=1    do not build/test the CUDA_DOUBLE version
#   CONTINUE_ON_SMOKE_FAILURE=1   run the full suites even if the smoke test fails
#   OLD_OXDNA        path to an original anm-oxdna CPU binary (scripts/build_old.sh): adds an
#                    "original CPU vs current CUDA" run of the suite
set -uo pipefail

if [ $# -lt 2 ]; then
    sed -n '2,30p' "$0"
    exit 2
fi
SRC="$(cd "$1" && pwd)"
OUT="$(mkdir -p "$2" && cd "$2" && pwd)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# under sbatch the script is copied to a spool directory: fall back to the submission directory
if [ ! -f "$REPO/run_tests.py" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
    [ -f "$REPO/run_tests.py" ] || REPO="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
[ -f "$REPO/run_tests.py" ] || { echo "cannot find run_tests.py (submit from the oxDNA_regression_tests directory)"; exit 2; }
PYTHON="${PYTHON:-python3}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 8)}}"
GPU_JOBS="${GPU_JOBS:-4}"
CMAKE_EXTRA="${CMAKE_EXTRA:--DCUDA_COMMON_ARCH=OFF}"
LOG="$OUT/cluster_cuda.log"
SUMMARY="$OUT/SUMMARY.md"
STATUS=0

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

if [ -n "${CLUSTER_MODULES:-}" ]; then
    type module > /dev/null 2>&1 || source /etc/profile.d/modules.sh 2> /dev/null || source /usr/share/lmod/lmod/init/bash 2> /dev/null
    # shellcheck disable=SC2086
    module load $CLUSTER_MODULES || { echo "module load $CLUSTER_MODULES failed"; exit 2; }
fi

: > "$LOG"
{
    echo "# oxDNA ANM/oxNA CUDA regression run"
    echo
    echo "- date: $(date)"
    echo "- host: $(hostname)"
    echo "- oxDNA source: $SRC ($(git -C "$SRC" rev-parse --abbrev-ref HEAD 2>/dev/null) $(git -C "$SRC" rev-parse --short HEAD 2>/dev/null))"
    echo "- test repo: $REPO ($(git -C "$REPO" rev-parse --short HEAD 2>/dev/null))"
    echo "- GPU: $(nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader 2>/dev/null | head -1)"
    echo "- nvcc: $(nvcc --version 2>/dev/null | tail -1)"
    echo "- cmake: $(cmake --version 2>/dev/null | head -1)"
    echo
} > "$SUMMARY"
"$PYTHON" -c "import numpy" || { echo "PYTHON=$PYTHON has no numpy" | tee -a "$SUMMARY"; exit 2; }

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
        log "BUILD FAILED: $name (see $dir/cmake.log, $dir/make.log)"
        echo "- **build $name FAILED** (see build/$name/make.log)" >> "$SUMMARY"
        tail -30 "$dir/make.log" >> "$LOG" 2>/dev/null
        return 1
    fi
}

suite() {  # name, run_tests.py args...
    local name="$1"; shift
    log "running $name: run_tests.py $*"
    "$PYTHON" "$REPO/run_tests.py" --jobs "$JOBS" --gpu-jobs "$GPU_JOBS" --out "$OUT/results/$name" "$@" > "$OUT/results_$name.log" 2>&1
    local rc=$?
    local res="$OUT/results/$name/summary.md"
    {
        echo "## $name (exit code $rc)"
        echo
        echo "\`run_tests.py $*\`"
        echo
        if [ -f "$res" ]; then
            grep -E '^totals:|unexpected failure' "$res"
            echo
            if grep -qE '\| (FAIL|ERROR|XPASS) \|' "$res"; then
                echo "| case | check | status | details |"
                echo "|---|---|---|---|"
                grep -E '\| (FAIL|ERROR|XPASS) \|' "$res"
                echo
            fi
        else
            echo "no summary written, see results_$name.log"
            tail -5 "$OUT/results_$name.log"
        fi
        echo
    } >> "$SUMMARY"
    log "$name finished with exit code $rc"
    return $rc
}

CPU="$OUT/build/cpu/bin/oxDNA"
MIXED="$OUT/build/cuda_mixed/bin/oxDNA"
DOUBLE="$OUT/build/cuda_double/bin/oxDNA"

# 1. builds needed for the smoke test
build cpu -DCUDA=OFF || exit 1
build cuda_mixed -DCUDA=ON || { STATUS=1; }

# 2. smoke test (energies + a short NVE run on a few systems, ~1-2 minutes): a broken build fails here
if [ -x "$MIXED" ]; then
    suite smoke_mixed --same-code --old "$CPU" --new "$MIXED" --backend CUDA --cuda-precision mixed --smoke
    rc=$?
    if [ $rc -ne 0 ]; then
        STATUS=1
        if [ "${CONTINUE_ON_SMOKE_FAILURE:-0}" != "1" ]; then
            log "smoke test failed: stopping (set CONTINUE_ON_SMOKE_FAILURE=1 to run the full suites anyway)"
        fi
    fi
fi

if [ $STATUS -eq 0 ] || [ "${CONTINUE_ON_SMOKE_FAILURE:-0}" = "1" ]; then
    # 3. CUDA double precision build (use_edge = false: the neighbour-matrix kernels)
    if [ "${SKIP_DOUBLE:-0}" != "1" ]; then
        build cuda_double -DCUDA=ON -DCUDA_DOUBLE=ON || STATUS=1
    fi

    # 4. full suites
    [ -x "$MIXED" ] && { suite cuda_mixed --same-code --old "$CPU" --new "$MIXED" --backend CUDA --cuda-precision mixed || STATUS=1; }
    if [ "${SKIP_DOUBLE:-0}" != "1" ] && [ -x "$DOUBLE" ]; then
        suite smoke_double --same-code --old "$CPU" --new "$DOUBLE" --backend CUDA --cuda-precision double --smoke || STATUS=1
        suite cuda_double --same-code --old "$CPU" --new "$DOUBLE" --backend CUDA --cuda-precision double || STATUS=1
    fi
    if [ -n "${OLD_OXDNA:-}" ] && [ -x "$MIXED" ]; then
        suite original_vs_cuda_mixed --old "$OLD_OXDNA" --new "$MIXED" --backend CUDA --cuda-precision mixed || STATUS=1
    fi

    # 5. the oxDNA repository's own CUDA NA tests (test/CUDA/NA, PR #230), if present
    if [ -f "$SRC/test/CUDA/NA/test_cuda_na.py" ]; then
        for prec in mixed double; do
            bin="$MIXED"; [ "$prec" = double ] && bin="$DOUBLE"
            [ -x "$bin" ] || continue
            log "test/CUDA/NA/test_cuda_na.py ($prec)"
            OXDNA_CUDA_PRECISION=$prec "$PYTHON" "$SRC/test/CUDA/NA/test_cuda_na.py" "$bin" -v > "$OUT/test_cuda_na_$prec.log" 2>&1
            rc=$?
            [ $rc -ne 0 ] && STATUS=1
            echo "## oxDNA test/CUDA/NA/test_cuda_na.py ($prec): exit code $rc" >> "$SUMMARY"
            echo >> "$SUMMARY"
            echo '```' >> "$SUMMARY"
            tail -4 "$OUT/test_cuda_na_$prec.log" >> "$SUMMARY"
            echo '```' >> "$SUMMARY"
            echo >> "$SUMMARY"
        done
    fi
fi

echo "**overall: $([ $STATUS -eq 0 ] && echo PASS || echo FAIL)**" >> "$SUMMARY"

# 6. tarball: summaries, results.json, inputs and logs of every run (trajectories left out to keep it small)
TAR="$OUT/oxdna_anm_cuda_results_$(date +%Y%m%d_%H%M).tar.gz"
(cd "$OUT" && find . -path ./build -prune -o -type f \( -name 'summary.md' -o -name 'results.json' -o -name 'input' \
    -o -name 'log.dat' -o -name 'stdout.txt' -o -name 'energy.dat' -o -name 'ext.dat' -o -name '*.log' -o -name 'SUMMARY.md' \) -print \
    | tar -czf "$TAR" -T - 2>/dev/null; tar -czf build_logs.tar.gz build/*/cmake.log build/*/make.log 2>/dev/null)
log "summary: $SUMMARY"
log "tarball: $TAR (build logs: $OUT/build_logs.tar.gz)"
cat "$SUMMARY"
exit $STATUS
