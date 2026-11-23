#!/usr/bin/env julia
# Run the 1-D reference of the porous shell validation over formation times and initial
# porosities (R_p 20 km, Nr 2000, dt 2 kyr, to the shell freeze or 20 Myr).
# Usage: julia --project=. scripts/run_validation_porous_shell_table.jl <out.csv> t_form:phi_init ...
# Example: ... table.csv 1.5:0.25 2.5:0.25 (t_form in Myr). A case takes 3 to 110 min.

include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
using .PorousShellValidation: solve_porous_sphere

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 ||
        error("Usage: run_validation_porous_shell_table.jl <out.csv> t:phi ...")
    cases = [parse.(Float64, split(c, ":")) for c in ARGS[2:end]]
    open(ARGS[1], "w") do io
        println(
            io,
            "t_form_myr,phi_init,T_peak,t_peak_myr,t_stop_myr,d_shell_m,d_phi04_m,eres,n_steps,wall_s",
        )
        for (tf, ph) in cases
            t = @elapsed r = solve_porous_sphere(;
                R=2.0e4,
                phi_init=ph,
                t_form_yr=tf * 1.0e6,
                t_end_yr=2.0e7,
                dt_yr=2.0e3,
                Nr=2000,
            )
            row = (
                tf,
                ph,
                round(r.T_peak; digits=1),
                round(r.t_peak_yr / 1.0e6; digits=3),
                round(r.t_stop_yr / 1.0e6; digits=3),
                round(r.d_shell; digits=1),
                round(r.d_phi04; digits=1),
                r.energy_residual,
                r.n_steps,
                round(t; digits=1),
            )
            println(io, join(row, ","))
            flush(io)
        end
    end
end
