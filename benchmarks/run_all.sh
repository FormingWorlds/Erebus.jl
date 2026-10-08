#!/usr/bin/env bash
# Execute all benchmark exporters and figure generators for Erebus.jl.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== Erebus.jl Benchmarks Suite Runner ==="
echo "Repository root: $REPO_ROOT"
echo "Benchmarks dir:  $SCRIPT_DIR"

mkdir -p "$REPO_ROOT/output_files"
mkdir -p "$REPO_ROOT/docs/src/assets"

echo ""
echo "--- Running Julia Benchmark Exporters ---"
julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/export_gravity_two_layer_benchmark.jl"
julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/export_hydrofracture_ramp_benchmark.jl"
julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/export_jeans_effusion_benchmark.jl"
julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/export_radiogenic_decay_benchmark.jl"
julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/export_thermal_slab_benchmark.jl"
julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/export_degassing_benchmark.jl"

echo ""
echo "--- Running Python Benchmark Figure Generators ---"
python3 "$SCRIPT_DIR/generate_accretion_benchmarks.py"
python3 "$SCRIPT_DIR/generate_core_formation_benchmark.py"
python3 "$SCRIPT_DIR/generate_core_geochemistry_benchmarks.py"
python3 "$SCRIPT_DIR/generate_coupled_atmosphere_benchmarks.py"
python3 "$SCRIPT_DIR/generate_degassing_benchmark.py"
python3 "$SCRIPT_DIR/generate_gravity_two_layer_benchmark.py"
python3 "$SCRIPT_DIR/generate_hydrofracture_ramp_benchmark.py"
python3 "$SCRIPT_DIR/generate_hydrothermal_convection_benchmarks.py"
python3 "$SCRIPT_DIR/generate_jeans_effusion_benchmark.py"
python3 "$SCRIPT_DIR/generate_lunar_growth_benchmarks.py"
python3 "$SCRIPT_DIR/generate_mineral_assemblage_benchmarks.py"
python3 "$SCRIPT_DIR/generate_multistage_accretion_benchmarks.py"
python3 "$SCRIPT_DIR/generate_radiogenic_decay_benchmark.py"
python3 "$SCRIPT_DIR/generate_reservoir_diagram.py"
python3 "$SCRIPT_DIR/generate_telescoping_benchmarks.py"
python3 "$SCRIPT_DIR/generate_thermal_slab_benchmark.py"
python3 "$SCRIPT_DIR/generate_volatile_mixture_benchmarks.py"
python3 "$SCRIPT_DIR/plot_diagnostics.py"
python3 "$SCRIPT_DIR/render_core_formation_movie.py"
python3 "$SCRIPT_DIR/render_magma_ocean_movie.py"

echo ""
echo "=== All benchmark exporters and generators completed successfully ==="
