#!/bin/bash
set -e
echo "Starting porosity parameter sweep (phi = 0.20, 0.35, 0.50)..."
julia --project=. -e 'using Erebus; cfg = load_config("configs/hydrothermal_reaction_sweep_phi20.toml"); Erebus.simulation_loop(cfg)' > /tmp/sweep_phi20.log 2>&1 &
PID20=$!
julia --project=. -e 'using Erebus; cfg = load_config("configs/hydrothermal_reaction_sweep_phi35.toml"); Erebus.simulation_loop(cfg)' > /tmp/sweep_phi35.log 2>&1 &
PID35=$!
julia --project=. -e 'using Erebus; cfg = load_config("configs/hydrothermal_reaction_sweep_phi50.toml"); Erebus.simulation_loop(cfg)' > /tmp/sweep_phi50.log 2>&1 &
PID50=$!
trap 'kill $PID20 $PID35 $PID50 2>/dev/null || true' EXIT

echo "Launched phi=0.20 (PID $PID20), phi=0.35 (PID $PID35), phi=0.50 (PID $PID50)"
wait $PID20
echo "phi=0.20 complete!"
wait $PID35
echo "phi=0.35 complete!"
wait $PID50
echo "phi=0.50 complete!"
echo "All sweep simulations completed successfully."
