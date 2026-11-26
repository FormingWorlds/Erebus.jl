#!/usr/bin/env julia
# Run the melting shell case (Erebus against the 1-D reference) on the grids and latent
# heats of RUNS and write one result row and one profile per run into a folder.
# Usage: julia --project=. scripts/run_validation_melting_shell.jl <out_dir> [tag ...]
# Run time on one core: about 2 min per run at Nx 65 and 7 min at Nx 129.

include(joinpath(@__DIR__, "..", "benchmarks", "validation_melting_shell.jl"))
using .MeltingShellValidation:
    run_melting_shell_case,
    reference_solution,
    L_MELT,
    L_MUTANT_FACTOR,
    NR_REFERENCE,
    PHI_INIT_CASE,
    T_FORM_YR_CASE,
    N_STEPS_CASE,
    DT_YR_CASE,
    SUITE_RESULT_HEADER

# (tag, Nx, latent heat of the Erebus run [J/kg]); every run compares with the 1-D
# reference at the nominal latent heat L_MELT.
const RUNS = (
    ("nx65", 65, L_MELT),
    ("nx129", 129, L_MELT),
    ("nx65_L0", 65, 0.0),
    ("nx65_L1.5", 65, L_MUTANT_FACTOR * L_MELT),
)

"""
    run_case(outdir, run, ref) -> String

Run one entry `(tag, Nx, L_melt)` of `RUNS` against the 1-D result `ref`, write
`result_<tag>.csv` (one row of `SUITE_RESULT_HEADER`) and `profile_<tag>.csv` (binned
peak temperatures of Erebus and of the 1-D reference) into `outdir`, and return the row.
"""
function run_case(outdir::String, run, ref)
    tag, nx, L = run
    t = @elapsed r = run_melting_shell_case(;
        Nx=nx,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        n_steps=N_STEPS_CASE,
        dt_yr=DT_YR_CASE,
        L_melt=L,
        Nr=NR_REFERENCE,
        reference=ref,
    )
    open(joinpath(outdir, "profile_$tag.csv"), "w") do io
        println(io, "r_m,tmax_K,tmax_ref_K")
        for i in eachindex(r.r)
            println(io, join((r.r[i], r.tmax[i], r.tmax_ref[i]), ","))
        end
    end
    row = join(
        (
            tag,
            nx,
            L,
            r.dx,
            r.T_centre,
            r.T_centre_ref,
            r.r_melt,
            r.r_melt_ref,
            r.tmax_rel_l2,
            r.energy_lse_residual,
            r.enthalpy_residual_ref,
            round(t; digits=1),
        ),
        ",",
    )
    write(joinpath(outdir, "result_$tag.csv"), row * "\n")
    return row
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: run_validation_melting_shell.jl <out_dir> [tag ...]")
    outdir = mkpath(ARGS[1])
    write(joinpath(outdir, "header.csv"), join(SUITE_RESULT_HEADER, ",") * "\n")
    ref = reference_solution()
    open(joinpath(outdir, "profile_1d.csv"), "w") do io
        println(io, "r_m,tmax_K")
        for i in eachindex(ref.r)
            println(io, ref.r[i], ",", ref.T_max[i])
        end
    end
    tags = length(ARGS) > 1 ? ARGS[2:end] : [r[1] for r in RUNS]
    for run in RUNS
        if run[1] in tags
            println(run_case(outdir, run, ref))
            flush(stdout)
            GC.gc()
        end
    end
end
