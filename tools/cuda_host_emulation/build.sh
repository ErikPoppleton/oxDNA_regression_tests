#!/usr/bin/env bash
# Checks the CUDA ANM code without a GPU.
# usage: build.sh <oxDNA checkout> <strict CPU build dir (cmake with CMAKE_CXX_COMPILER_LAUNCHER adding -fno-fast-math -ffp-contract=off)>
#  1. parses the CUDA sources with clang (host and device side, float and double) against the stub CUDA API in mock/
#  2. links the host objects of CUDAANMInteraction.cu and CUDADRHInteraction.cu together (duplicate symbols)
#  3. builds test_anm_device (device functions vs ANMModel) and test_anm_kernels (whole kernels vs the CPU interaction)
# Run: ./test_anm_kernels DNACT <topology> <conf> <parfile> [kb kt]
set -e
SRC=$1/src; BUILD=$2; HERE=$(cd "$(dirname "$0")" && pwd)
CLANG=${CLANG:-/opt/homebrew/opt/llvm/bin/clang++}
for f in CUDA/Interactions/CUDAANMInteraction.cu CUDA/Interactions/CUDADRHInteraction.cu CUDA/Interactions/CUDAInteractionFactory.cu; do
	for def in "" "-DCUDA_DOUBLE_PRECISION"; do
		for mode in --cuda-host-only --cuda-device-only; do
			(cd $SRC && $CLANG -x cuda $mode $def -nocudainc -nocudalib --cuda-gpu-arch=sm_70 -include $HERE/mock/cuda_runtime_api.h -I$HERE/mock -Iextern -std=c++14 -fsyntax-only -Wno-everything $f) && echo "OK $f $def $mode"
		done
	done
done
# 2. link the host objects of the .cu files that share device code into one library, as the real build does:
#    kernels or device functions defined in a header included by two of them are then defined twice
OBJ=$(mktemp -d)
for def in "" "-DCUDA_DOUBLE_PRECISION"; do
	for f in CUDAANMInteraction CUDADRHInteraction; do
		(cd $SRC && $CLANG -x cuda --cuda-host-only $def -nocudainc -nocudalib --cuda-gpu-arch=sm_70 -include $HERE/mock/cuda_runtime_api.h -I$HERE/mock -Iextern -std=c++14 -Wno-everything -c CUDA/Interactions/$f.cu -o $OBJ/$f.o)
	done
	if [ "$(uname)" = Darwin ]; then UNDEF="-Wl,-undefined,dynamic_lookup"; else UNDEF=""; fi
	$CLANG -shared $UNDEF $OBJ/CUDAANMInteraction.o $OBJ/CUDADRHInteraction.o -o $OBJ/liblinkcheck.so && echo "OK link CUDAANMInteraction + CUDADRHInteraction $def"
done
rm -rf $OBJ
python3 - "$SRC" "$HERE" <<'PY'
import sys
src, here = sys.argv[1:]
s = open(src + "/CUDA/Interactions/CUDARNAInteraction.cu").read()
i = s.index("void copy_Model_to_CUDAModel(Model& model_from, CUDAModel& model_to) {")
j = s.index("\n}\n", i) + 3
body = s[i:j].replace("void copy_Model_to_CUDAModel(Model& model_from, CUDAModel& model_to)", "static void copy_impl(Model& model_from, cuda_rna::CUDAModel& model_to)")
body += "void copy_Model_to_CUDAModel(Model &model_from, CUDAModel &model_to) { copy_impl(model_from, reinterpret_cast<cuda_rna::CUDAModel &>(model_to)); }\n"
open(here + "/copy_model.inc", "w").write(body)
PY
TE=$(find $BUILD -name "libtinyexpr*" | head -1)
for t in test_anm_device test_anm_kernels; do
	(cd $SRC && clang++ -std=c++14 -DCUDA_DOUBLE_PRECISION -DHOST_EMULATION -fno-fast-math -ffp-contract=off -I$HERE/mock -I$HERE -Iextern -I. -O2 -Wno-everything $HERE/$t.cpp $HERE/stubs.cpp $BUILD/src/liboxdna_common.a $TE -o $HERE/$t)
done
echo "built $HERE/test_anm_device and $HERE/test_anm_kernels"
