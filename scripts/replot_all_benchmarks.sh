#!/usr/bin/env bash
set -euo pipefail

# replot_all_benchmarks.sh
# Master script to regenerate all documentation benchmark and verification figures for Erebus.jl.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "=== Regenerating All Documentation Benchmark Figures ==="

# Select Python 3 interpreter
if command -v python3 >/dev/null 2>&1; then
    PYTHON="python3"
elif [ -x "/Users/timlichtenberg/miniforge3/bin/python3" ]; then
    PYTHON="/Users/timlichtenberg/miniforge3/bin/python3"
else
    echo "Error: python3 interpreter not found." >&2
    exit 1
fi

echo "Using Python: $(${PYTHON} --version) at $(which ${PYTHON})"

# 1. Poroelastic constitutive limits verification
echo "[1/10] Poroelastic constitutive limits..."
"${PYTHON}" "${SCRIPT_DIR}/generate_poroelastic_benchmark.py"

# 2. 1D Terzaghi analytical consolidation benchmark
echo "[2/10] Terzaghi consolidation benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_terzaghi_benchmark.py"

# 3. Comprehensive H-C-N-S volatile solubility benchmark
echo "[3/10] H-C-N-S volatile solubility benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_hcns_solubility_benchmark.py"

# 4. Volatile solubility & redox-dependent nitrogen benchmark
echo "[4/10] Volatile solubility benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_volatile_solubility_benchmark.py"

# 5. Cold lid hydrofracture and episodic venting benchmark
echo "[5/10] Cold surface venting benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_cold_venting_benchmark.py"

# 6. Hydrofracture regime and cryogenic sealing benchmark
echo "[6/10] Hydrofracture venting benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_hydrofracture_venting_benchmark.py"

# 7. Atmospheric Jeans kinetic escape benchmark
echo "[7/10] Jeans escape benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_jeans_escape_benchmark.py"

# 8. Darcy buoyancy verification benchmark
echo "[8/10] Darcy buoyancy verification..."
"${PYTHON}" "${SCRIPT_DIR}/generate_buoyancy_benchmark.py"

# 9. Fluid viscosity temperature dependence benchmark
echo "[9/10] Fluid viscosity benchmark..."
"${PYTHON}" "${SCRIPT_DIR}/generate_viscosity_benchmark.py"

# 10. Hydrofracture tensile permeability verification
echo "[10/10] Hydrofracture verification..."
"${PYTHON}" "${SCRIPT_DIR}/generate_hydrofracture_benchmark.py"

echo "=== All benchmark figures successfully regenerated in docs/src/assets/ ==="
