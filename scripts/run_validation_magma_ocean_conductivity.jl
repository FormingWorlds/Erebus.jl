#!/usr/bin/env julia
# Run the magma-ocean conductivity cases against the 1-D reference and write CSV rows.
# Usage: julia --project=. scripts/run_validation_magma_ocean_conductivity.jl <out_dir> [tag ...]

include(joinpath(@__DIR__, "..", "benchmarks", "validation_magma_ocean_conductivity.jl"))
using .MagmaOceanConductivityValidation:
    PSV,
    MSV,
    run_magma_ocean_case,
    reference_solution,
    reference_soft_turb,
    TURB_EXP_REF,
    TURB_EXP_LOW_MUTANT,
    TURB_EXP_HIGH_MUTANT,
    NR_REFERENCE,
    PHI_INIT_CASE,
    T_FORM_YR_CASE,
    N_STEPS_CASE,
    DT_YR_CASE,
    SUITE_RESULT_HEADER

# (tag, Nx, soft_turbulence, turb_exponent); every run compares against the nominal 1-D model.
const RUNS = (
    ("nx65", 65, true, TURB_EXP_REF),
    ("nx129", 129, true, TURB_EXP_REF),
    ("nx65_off", 65, false, TURB_EXP_REF),
    ("nx65_exp02", 65, true, TURB_EXP_LOW_MUTANT),
    ("nx65_exp05", 65, true, TURB_EXP_HIGH_MUTANT),
)

"""
    run_case(outdir, run, ref) -> String

Run one entry `(tag, Nx, soft_turbulence, turb_exponent)` of `RUNS` against `ref`, write
`result_<tag>.csv` and `profile_<tag>.csv` into `outdir`, and return the CSV row.
"""
function run_case(outdir::String, run, ref)
    tag, nx, st_flag, bexp = run
    t = @elapsed r = run_magma_ocean_case(;
        Nx=nx,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        n_steps=N_STEPS_CASE,
        dt_yr=DT_YR_CASE,
        soft_turbulence=st_flag,
        turb_exponent=bexp,
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
            st_flag,
            bexp,
            r.dx,
            r.T_centre,
            r.T_centre_ref,
            r.plateau_drop,
            r.plateau_drop_ref,
            r.r_mush,
            r.r_mush_ref,
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
    isempty(ARGS) &&
        error("Usage: run_validation_magma_ocean_conductivity.jl <out_dir> [tag ...]")
    outdir = mkpath(ARGS[1])
    hdr_tmp = joinpath(outdir, ".header_$(getpid()).tmp")
    write(hdr_tmp, join(SUITE_RESULT_HEADER, ",") * "\n")
    mv(hdr_tmp, joinpath(outdir, "header.csv"); force=true)
    ref = reference_solution()
    ref_off = PSV.solve_porous_sphere(;
        R=PSV.R_PLANET,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        t_end_yr=N_STEPS_CASE * DT_YR_CASE,
        dt_yr=DT_YR_CASE,
        Nr=NR_REFERENCE,
        stop_cold=false,
        melt=MSV.reference_melt(),
        soft_turb=reference_soft_turb(; enabled=false),
    )
    p1d_tmp = joinpath(outdir, ".profile_1d_$(getpid()).tmp")
    open(p1d_tmp, "w") do io
        println(io, "r_m,tmax_K,tmax_off_K")
        for i in eachindex(ref.r)
            println(io, join((ref.r[i], ref.T_max[i], ref_off.T_max[i]), ","))
        end
    end
    mv(p1d_tmp, joinpath(outdir, "profile_1d.csv"); force=true)
    tags = length(ARGS) > 1 ? ARGS[2:end] : [r[1] for r in RUNS]
    for run in RUNS
        if run[1] in tags
            println(run_case(outdir, run, ref))
            flush(stdout)
            GC.gc()
        end
    end
end
