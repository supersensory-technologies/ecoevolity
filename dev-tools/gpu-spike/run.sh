#!/usr/bin/env bash
# Build and run the GPU spike benchmark (expmv_bench.cu): phycoeval's per-branch
# matrix-exponential step in FP64 on an NVIDIA GPU, checked against
# ecoevolity's own CPU code. Needs nvcc and a GPU; ARCH defaults to sm_120.
#
#   dev-tools/gpu-spike/run.sh [patterns]
#
# On a shared GPU host, run it under the host's GPU lock, e.g.
#   flock /run/lock/sst-gpu.lock dev-tools/gpu-spike/run.sh
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src="$here/../../src/ecoevolity"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The 12 CF poles and residues, copied out of MatrixExponentiator.
awk '
/const std::array<double, 12> (ci_real_|ci_imag_|zi_real_|zi_imag_)/ {
    match($0, /(ci_real_|ci_imag_|zi_real_|zi_imag_)/); name = substr($0, RSTART, RLENGTH); f = 1; vals = ""; next
}
f && /}/ {
    map["ci_real_"] = "ECO_CR"; map["ci_imag_"] = "ECO_CI"; map["zi_real_"] = "ECO_ZR"; map["zi_imag_"] = "ECO_ZI"
    print "extern const double " map[name] "[12] = {" vals "};"; f = 0; next
}
f { gsub(/[ \t]/, ""); if ($0 != "") vals = vals $0 }
' "$src/matrix.hpp" > "$work/consts.cpp"

# --fmad=false keeps the GPU's arithmetic identical to the CPU's.
nvcc -O2 -std=c++17 -arch="${ARCH:-sm_120}" --fmad=false -I"$src" \
    -o "$work/expmv_bench" "$here/expmv_bench.cu" "$work/consts.cpp" "$src/error.cpp"
"$work/expmv_bench" "${1:-477}"
