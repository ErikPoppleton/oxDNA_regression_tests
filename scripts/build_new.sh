#!/usr/bin/env bash
# Build the current oxDNA at a given commit from a detached git worktree (the checkout itself is
# never touched, so this is safe while someone else is working in it).
#
#   scripts/build_new.sh [COMMIT] [WORKTREE_DIR] [extra cmake args...]
#
#   COMMIT        commit/branch/tag to build      (default: $NEW_COMMIT or f802bea4)
#   WORKTREE_DIR  where to put the worktree       (default: <repo>/build/new_<commit>)
#   NEW_REPO      oxDNA repository (env)          (default: ~/software/oxDNA)
#   e.g. CUDA: scripts/build_new.sh f802bea4 build/new_cuda -DCUDA=ON
#
# CPU builds are double precision by default (DOUBLE=ON); CUDA builds use mixed precision.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
NEW_REPO="${NEW_REPO:-$HOME/software/oxDNA}"
COMMIT="${1:-${NEW_COMMIT:-f802bea4}}"
WT="${2:-$REPO/build/new_$COMMIT}"
shift $(( $# > 2 ? 2 : $# )) || true
JOBS="${JOBS:-8}"

if [ ! -d "$WT" ]; then
    git -C "$NEW_REPO" worktree add --detach "$WT" "$COMMIT"
else
    echo "re-using existing worktree $WT (at $(git -C "$WT" rev-parse --short HEAD))"
fi
mkdir -p "$WT/build"
cmake -S "$WT" -B "$WT/build" "$@"
make -C "$WT/build" -j"$JOBS" oxDNA DNAnalysis
echo "new binary: $WT/build/bin/oxDNA"
