#!/usr/bin/env julia
# Run the full-evolution cases of the porous shell validation (Erebus against the 1-D
# reference) and write one result row and one profile per case into a folder.
# Usage: julia --project=. scripts/run_validation_porous_shell_suite.jl <out_dir> [tag ...]
# A case takes 13 to 45 min at Nx 129 on one core; run at most a few at once.

include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
using .PorousShellValidation:
    run_porous_shell_case, NR_REFERENCE, SUITE_CASES, SUITE_DT_YR, SUITE_RESULT_HEADER

# The case list is `SUITE_CASES` of the validation module: (tag, Nx, t_form [Myr],
# phi_init, steps of `SUITE_DT_YR` to the shell freeze, air conductivity).
const RESULT_HEADER = join(SUITE_RESULT_HEADER, ",")

"""
    run_suite_case(outdir, case) -> String

Run one entry of `SUITE_CASES`, write `result_<tag>.csv` (one row) and `prof_<tag>.csv`
(binned profiles) into `outdir`, and return the result row.
"""
function run_suite_case(outdir::String, case)
    tag, Nx, tf, ph, n, kair = case
    t = @elapsed r = run_porous_shell_case(;
        Nx=Nx,
        phi_init=ph,
        t_form_yr=tf * 1.0e6,
        n_steps=n,
        dt_yr=SUITE_DT_YR,
        k_air=kair,
        Nr=NR_REFERENCE,
    )
    open(joinpath(outdir, "prof_$tag.csv"), "w") do io
        println(io, "r_m,tmax_K,tmax_ref_K,phi,phi_ref")
        for i in eachindex(r.r)
            println(
                io, join((r.r[i], r.tmax[i], r.tmax_ref[i], r.phi[i], r.phi_ref[i]), ",")
            )
        end
    end
    row = join(
        (
            tag,
            Nx,
            tf,
            ph,
            n,
            something(kair, -1.0),
            r.dx,
            r.d_shell,
            r.d_shell_q25,
            r.d_shell_q75,
            r.d_shell_ref,
            round(r.d_shell_error_cells; digits=3),
            round(r.tmax_rel_l2; digits=4),
            round(r.phi_rms; digits=4),
            round(r.T_peak; digits=2),
            round(r.T_peak_ref; digits=2),
            r.energy_lse_residual,
            round(t; digits=1),
        ),
        ",",
    )
    open(joinpath(outdir, "result_$tag.csv"), "w") do io
        return println(io, row)
    end
    return row
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("Usage: run_validation_porous_shell_suite.jl <out_dir> [tag ...]")
    outdir = mkpath(ARGS[1])
    write(joinpath(outdir, "header.csv"), RESULT_HEADER * "\n")
    tags = length(ARGS) > 1 ? ARGS[2:end] : [c[1] for c in SUITE_CASES]
    for case in SUITE_CASES
        case[1] in tags && println(run_suite_case(outdir, case))
    end
end
