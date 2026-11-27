using Test
using TOML
using StaticArrays
using Erebus

isdefined(@__MODULE__, :MagmaOceanConductivityValidation) || include(
    joinpath(@__DIR__, "..", "benchmarks", "validation_magma_ocean_conductivity.jl")
)
isdefined(@__MODULE__, :run_validation_cli) ||
    include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
isdefined(@__MODULE__, :MagmaOceanConductivityExport) || include(
    joinpath(@__DIR__, "..", "tools", "export_validation_magma_ocean_conductivity.jl")
)

const MOSV = MagmaOceanConductivityValidation
const MO_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "magma_ocean_conductivity.toml"
)
const MO_SUITE_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "magma_ocean_conductivity_suite.toml"
)

"""
Write the record dictionary `dict` as TOML to `path` and return `path`.
"""
function _mo_write_toml(path::String, dict::AbstractDict)
    open(io -> TOML.print(io, dict), path, "w")
    return path
end

"""
Call `f()` with standard output and error redirected to temporary files.
"""
function _mo_capture(f)
    return mktempdir() do dir
        out_path, err_path = joinpath(dir, "out.txt"), joinpath(dir, "err.txt")
        value = open(out_path, "w") do out
            open(err_path, "w") do err
                return redirect_stdio(f; stdout=out, stderr=err)
            end
        end
        return (; value, out=read(out_path, String), err=read(err_path, String))
    end
end

"""
Initialize a simulation state and close its telemetry stream.
"""
function _mo_state(cfg)
    st, coords, ws, tel, _ = Erebus.init_simulation(cfg; output_path=mktempdir())
    tel !== nothing && close(tel)
    return st, coords, ws
end

"""
Set rock markers to temperature `T` and dry porosity `phi_d`, then interpolate to the grid.
"""
function _mo_marker_step!(st, coords, ws, cfg, T, phi_d)
    mk = st.markers
    rock = findall(<(3), mk.core.tm)
    mk.core.tkm[rock] .= T
    mk.groups.porosity.phi_d[rock] .= phi_d
    Erebus.radiogenic_heating!(st, coords, cfg)
    Erebus.interpolate_markers_to_grid!(
        st,
        coords,
        cfg;
        thread_buffers=ws.thread_buffers,
        interp_arrays=ws.interp_arrays,
        rock_sums=ws.rock_sums,
    )
    return rock
end

@testset "Magma Ocean Effective Conductivity Validation" begin
    @testset "Marker conductivity, grid faces and subgrid relaxation across melt regimes" begin
        kw = (; Nx=17, phi_init=0.25, t_form_yr=0.85e6, n_steps=1, dt_yr=1.0e4)
        cfg_on = MOSV.make_magma_ocean_config(; kw..., soft_turbulence=true)
        cfg_off = MOSV.make_magma_ocean_config(; kw..., soft_turbulence=false)
        @test cfg_on.melting.active && cfg_on.melting.soft_turbulence
        @test cfg_off.melting.active && !cfg_off.melting.soft_turbulence

        # Below F_turb_start = 0.30 (T = 1500 K, F = 0.25): soft turbulence is inactive
        st_on, coords, ws_on = _mo_state(cfg_on)
        st_off, _, ws_off = _mo_state(cfg_off)
        rock_on = _mo_marker_step!(st_on, coords, ws_on, cfg_on, 1500.0, 0.25)
        rock_off = _mo_marker_step!(st_off, coords, ws_off, cfg_off, 1500.0, 0.25)
        k_por_025 = MOSV.PSV.k_porous_ref(3.0, 0.25)
        @test all(isapprox.(st_on.markers.core.ktotalm[rock_on], k_por_025; rtol=1e-13))
        @test isapprox(
            st_on.markers.core.ktotalm[rock_on],
            st_off.markers.core.ktotalm[rock_off];
            rtol=1.0e-15,
            atol=0.0,
        )

        # Mush midpoint F = 0.40 (T = 1560 K, phi_d = 0): pins mush eta_fluid blend
        _mo_marker_step!(st_on, coords, ws_on, cfg_on, 1560.0, 0.0)
        eta_040 = 1.0e19 * exp(-28.0 * 0.40)
        eta_fl_040 = sqrt(eta_040 * 100.0)
        k_turb_040 = 3.0 * (eta_040 / eta_fl_040)^(1.0 / 3.0)
        k_exp_040 = sqrt(3.0 * k_turb_040)
        k_no_mush_blend = sqrt(3.0 * (3.0 * (eta_040 / 100.0)^(1.0 / 3.0)))
        k_laminar_040 = sqrt(3.0 * (3.0 * sqrt(eta_040 / eta_fl_040)))
        k_m040 = st_on.markers.core.ktotalm[rock_on[1]]
        @test isapprox(k_m040, k_exp_040; rtol=1e-12)
        @test isapprox(k_m040, 30.792603; atol=1e-5)
        @test 10.0 < k_m040 < 60.0
        @test 250.0 < k_no_mush_blend - k_m040 < 320.0
        @test 50.0 < k_laminar_040 - k_m040 < 85.0

        # Disaggregated suspension below etamin clamp F = 0.45 (T = 1580 K, phi_d = 0)
        _mo_marker_step!(st_on, coords, ws_on, cfg_on, 1580.0, 0.0)
        eta_045 = exp10(log10(eta_040) - (0.05 / 0.60) * (log10(eta_040) - 1.0))
        eta_fl_045 = exp(0.25 * log(eta_045) + 0.75 * log(100.0))
        k_turb_045 = 3.0 * (eta_045 / eta_fl_045)^(1.0 / 3.0)
        w_045 = 3.0 * 0.75^2 - 2.0 * 0.75^3
        k_exp_045 = exp10((1.0 - w_045) * log10(3.0) + w_045 * log10(k_turb_045))
        @test 5.0e12 < eta_045 < 2.0e13
        @test isapprox(st_on.markers.core.ktotalm[rock_on[1]], k_exp_045; rtol=1e-12)
        ic = (cfg_on.grid.Nx + 1) ÷ 2
        @test isapprox(st_on.grids.KX[ic, ic], k_exp_045; rtol=1e-12)
        @test isapprox(ws_on.rock_sums.KXu[ic, ic], k_exp_045; rtol=1e-12)

        # Suspension regime F = 0.625 (T = 1650 K): (etamin / eta_fluid_silicate)^(1/3)
        _mo_marker_step!(st_on, coords, ws_on, cfg_on, 1650.0, 0.0)
        k_susp = 3.0 * (1.0e12 / 100.0)^(1.0 / 3.0)
        @test isapprox(st_on.markers.core.ktotalm[rock_on[1]], k_susp; rtol=1e-12)
        @test isapprox(k_susp, 6463.304070; atol=1e-5)
        @test isapprox(st_on.grids.KX[ic, ic], k_susp; rtol=1e-12)
        @test isapprox(st_on.grids.KY[ic, ic], k_susp; rtol=1e-12)
        @test isapprox(ws_on.rock_sums.KXu[ic, ic], k_susp; rtol=1e-12)
        @test isapprox(ws_on.rock_sums.KYu[ic, ic], k_susp; rtol=1e-12)

        # Porosity factor on henke2012 in the suspension regime (phi_d = 0.25)
        _mo_marker_step!(st_on, coords, ws_on, cfg_on, 1650.0, 0.25)
        @test isapprox(
            st_on.markers.core.ktotalm[rock_on[1]],
            k_por_025 * (1.0e12 / 100.0)^(1.0 / 3.0);
            rtol=1e-12,
        )

        # Suspension parameter scaling: etamin = 1e13, eta_fluid_silicate = 1000, etasolidmm
        cfg_etamin = SimulationConfig(;
            [
                f => (
                    if f === :solver
                        SolverConfig(; seed=42, hydromech_solver=:direct, etamin=1.0e13)
                    else
                        getfield(cfg_on, f)
                    end
                ) for f in fieldnames(SimulationConfig)
            ]...,
        )
        st_em, c_em, ws_em = _mo_state(cfg_etamin)
        r_em = _mo_marker_step!(st_em, c_em, ws_em, cfg_etamin, 1650.0, 0.0)
        @test isapprox(st_em.markers.core.ktotalm[r_em[1]], k_susp * cbrt(10.0); rtol=1e-12)

        cfg_efl = MOSV.make_magma_ocean_config(; kw..., eta_fluid_silicate=1000.0)
        st_ef, c_ef, ws_ef = _mo_state(cfg_efl)
        r_ef = _mo_marker_step!(st_ef, c_ef, ws_ef, cfg_efl, 1650.0, 0.0)
        @test isapprox(st_ef.markers.core.ktotalm[r_ef[1]], k_susp / cbrt(10.0); rtol=1e-12)
        @test isapprox(st_ef.markers.core.ktotalm[r_ef[1]], 3000.0; rtol=1e-12)

        # Above tmsolidphase rock reads etasolidmm, not etasolidm
        mat_mm = MaterialConfig(;
            ksolidm=cfg_on.materials.ksolidm,
            rhocpsolidm=cfg_on.materials.rhocpsolidm,
            rhocpfluidm=cfg_on.materials.rhocpfluidm,
            alphasolidm=cfg_on.materials.alphasolidm,
            alphafluidm=cfg_on.materials.alphafluidm,
            tkm0=cfg_on.materials.tkm0,
            XWsolidm_init=cfg_on.materials.XWsolidm_init,
            etasolidm=SVector{3,Float64}([1.0e19, 1.0e19, 1.0e18]),
            etasolidmm=SVector{3,Float64}([1.0e15, 1.0e15, 1.0e18]),
        )
        cfg_mm = SimulationConfig(;
            [
                f => (f === :materials ? mat_mm : getfield(cfg_on, f)) for
                f in fieldnames(SimulationConfig)
            ]...,
        )
        st_mm, c_mm, ws_mm = _mo_state(cfg_mm)
        r_mm = _mo_marker_step!(st_mm, c_mm, ws_mm, cfg_mm, 1560.0, 0.0)
        k_mm_exp = sqrt(3.0 * 3.0 * (1.0e12 / sqrt(1.0e12 * 100.0))^(1.0 / 3.0))
        @test isapprox(st_mm.markers.core.ktotalm[r_mm[1]], k_mm_exp; rtol=1e-12)
        @test 18.0 < st_mm.markers.core.ktotalm[r_mm[1]] < 23.0

        # Subgrid relaxation uses ktotalm rather than k_porous
        _mo_marker_step!(st_on, coords, ws_on, cfg_on, 1650.0, 0.0)
        xc_on, yc_on = coords.xsize / 2.0, coords.ysize / 2.0
        m_id = rock_on[argmin(
            hypot.(
                st_on.markers.core.xm[rock_on] .- xc_on,
                st_on.markers.core.ym[rock_on] .- yc_on,
            ),
        )]
        st_on.markers.core.tkm[m_id] = 1660.0
        dt_sub = 1.0e9
        DT = zeros(coords.Ny1, coords.Nx1)
        TKSUM = zeros(coords.Ny1, coords.Nx1)
        RHOCPSUM = zeros(coords.Ny1, coords.Nx1)
        Erebus.apply_subgrid_temperature_diffusion!(
            st_on.markers.core.xm,
            st_on.markers.core.ym,
            st_on.markers.core.tm,
            st_on.markers.core.tkm,
            st_on.markers.core.phim,
            st_on.grids.tk1,
            DT,
            TKSUM,
            RHOCPSUM,
            dt_sub,
            length(st_on.markers),
            Erebus.marker_property_mode;
            coords=coords,
            dsubgridt=1.0,
            phi_dry=st_on.markers.groups.porosity.phi_d,
            rhocptotalm=st_on.markers.core.rhocptotalm,
            ktotalm=st_on.markers.core.ktotalm,
            fixed=st_on.grids.fixed_tk,
        )
        rhocp_m = st_on.markers.core.rhocptotalm[m_id]
        decay_turb = exp(-k_susp * dt_sub / rhocp_m * (4.0 / coords.dx^2))
        decay_cond = exp(-3.0 * dt_sub / rhocp_m * (4.0 / coords.dx^2))
        @test isapprox(st_on.markers.core.tkm[m_id] - 1650.0, 10.0 * decay_turb; rtol=1e-11)
        @test 2.0 < (10.0 * decay_cond) - (st_on.markers.core.tkm[m_id] - 1650.0) < 10.0
    end

    @testset "Explicit spherical metric step on an unrelaxed molten core jump" begin
        res = MOSV.hot_core_metric_step_check()
        @test isapprox(res.k_core_max, 6463.30; atol=0.05)
        @test isapprox(res.min_dT_metric, -1053.38; atol=0.1)
        @test isapprox(res.T_centre_2d_sph, 1489.57; atol=0.1)
        @test isapprox(res.T_centre_2d_no_metric, 1738.54; atol=0.1)
        @test isapprox(res.T_centre_1d, 1769.53; atol=0.1)
        @test isapprox(res.T_centre_2d_sph - res.T_centre_1d, -279.96; atol=0.1)
        @test -1200.0 < res.min_dT_metric < -1000.0
    end

    @testset "1-D soft-turbulence reference law and bitwise fallback" begin
        kw = (;
            R=2.0e4,
            phi_init=0.25,
            t_form_yr=1.15e6,
            t_end_yr=3.5e6,
            dt_yr=2.5e4,
            Nr=300,
            stop_cold=false,
            melt=MOSV.MSV.reference_melt(),
        )
        r_none = MOSV.PSV.solve_porous_sphere(; kw...)
        r_off = MOSV.PSV.solve_porous_sphere(;
            kw..., soft_turb=MOSV.reference_soft_turb(; enabled=false)
        )
        @test isapprox(r_none.T, r_off.T; rtol=1.0e-15, atol=0.0)
        @test isapprox(r_none.T_max, r_off.T_max; rtol=1.0e-15, atol=0.0)
        @test isapprox(r_none.phi, r_off.phi; rtol=1.0e-15, atol=0.0)

        # Reference k_soft_turb_ref matches compute_marker_properties! across regimes
        cfg = MOSV.make_case_config(; Nx=17)
        st, coords, ws = _mo_state(cfg)
        for pd in (0.0, 0.25),
            T_test in (1350.0, 1500.0, 1520.0, 1560.0, 1580.0, 1600.0, 1650.0, 1850.0)

            rock = _mo_marker_step!(st, coords, ws, cfg, T_test, pd)
            k_lat = MOSV.PSV.k_porous_ref(MOSV.PSV.K_SOL, pd)
            k_ref = MOSV.PSV.k_soft_turb_ref(
                k_lat, T_test; soft_turb=MOSV.reference_soft_turb()
            )
            @test isapprox(st.markers.core.ktotalm[rock[1]], k_ref; rtol=1e-13)
        end

        # Monotonicity in temperature and domain validation
        st_nt = MOSV.reference_soft_turb()
        ks = [MOSV.PSV.k_soft_turb_ref(3.0, T; soft_turb=st_nt) for T in 1400.0:25.0:1850.0]
        @test all(diff(ks) .>= 0.0)
        @test ks[end] > 2000.0 * ks[1]
        @test_throws DomainError MOSV.PSV.k_soft_turb_ref(0.0, 1600.0; soft_turb=st_nt)
        @test_throws DomainError MOSV.PSV.k_soft_turb_ref(-1.0, 1600.0; soft_turb=st_nt)
        @test_throws DomainError MOSV.PSV.k_soft_turb_ref(
            3.0, 1600.0; soft_turb=MOSV.reference_soft_turb(; turb_exponent=0.0)
        )
        @test_throws DomainError MOSV.PSV.k_soft_turb_ref(
            3.0, 300.0; soft_turb=MOSV.reference_soft_turb(; dT_min=0.0)
        )
        @test_throws DomainError MOSV.PSV.k_soft_turb_ref(
            3.0, 1600.0; soft_turb=MOSV.reference_soft_turb(; T_sol=1800.0, T_liq=1800.0)
        )
        @test_throws ArgumentError MOSV.PSV.solve_porous_sphere(;
            kw..., soft_turb=MOSV.reference_soft_turb(; face_mode=:harmonc)
        )
        @test_throws ArgumentError MOSV.PSV.solve_porous_sphere(;
            kw..., soft_turb=MOSV.reference_soft_turb(; T_sol=1450.0)
        )
        @test isapprox(
            MOSV.PSV.sinter_step_ref(1.0e-7, 0.25, 10.0, 1600.0, 1.0e11),
            1.0e-7;
            rtol=1.0e-15,
            atol=0.0,
        )
        @test_throws ArgumentError MOSV.make_magma_ocean_config(;
            Nx=17, phi_init=0.25, t_form_yr=0.85e6, n_steps=1, dt_yr=1.0e4, T_sol=1900.0
        )
    end

    @testset "Reference literal guards and canonical config SHA-256" begin
        ok, errs = MOSV.check_reference_literals_match_config()
        @test ok
        @test isempty(errs)
        @test isapprox(MOSV.T_MUSH, 1520.0; rtol=1e-14)

        # Existing porous_shell and melting_shell hashes stay untouched
        ps_raw = TOML.parsefile(
            joinpath(@__DIR__, "..", "benchmarks", "records", "porous_shell.toml")
        )
        ms_raw = TOML.parsefile(
            joinpath(@__DIR__, "..", "benchmarks", "records", "melting_shell.toml")
        )
        @test MOSV.PSV.canonical_config_sha256(MOSV.PSV.make_case_config()) ==
            String(ps_raw["metadata"]["config_sha256"])
        @test MOSV.MSV.canonical_config_sha256(MOSV.MSV.make_case_config()) ==
            String(ms_raw["metadata"]["config_sha256"])

        base = MOSV.make_case_config()
        sha = MOSV.canonical_config_sha256(base)
        @test length(sha) == 64
        @test MOSV.suite_config_sha256() != sha

        # Perturbing each rheology or soft-turbulence config field trips the guard and hash
        function with_sec(sec_name::Symbol, sec_val)
            return SimulationConfig(;
                [
                    f => (f === sec_name ? sec_val : getfield(base, f)) for
                    f in fieldnames(SimulationConfig)
                ]...,
            )
        end
        m0 = base.melting
        melt_mutants = (
            (
                "ALPHA_ETA_REF",
                MeltingConfig(;
                    [
                        f => (f === :alpha_eta ? 25.0 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "PHI_CRIT_REF",
                MeltingConfig(;
                    [
                        f => (f === :phi_crit ? 0.35 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "ETA_MELT_REF",
                MeltingConfig(;
                    [
                        f => (f === :eta_melt ? 20.0 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "ETA_FLUID_SIL_REF",
                MeltingConfig(;
                    [
                        f => (f === :eta_fluid_silicate ? 200.0 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "F_TURB_START_REF",
                MeltingConfig(;
                    [
                        f => (f === :F_turb_start ? 0.25 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "F_TURB_END_REF",
                MeltingConfig(;
                    [
                        f => (f === :F_turb_end ? 0.55 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "TURB_EXP_REF",
                MeltingConfig(;
                    [
                        f => (f === :turb_exponent ? 0.5 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "DT_TURB_MIN_REF",
                MeltingConfig(;
                    [
                        f => (f === :dT_turb_min ? 20.0 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "T_SURF_REF",
                MeltingConfig(;
                    [
                        f => (f === :T_surface_ref ? 310.0 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "K_TURB_FLOOR_REF",
                MeltingConfig(;
                    [
                        f => (f === :k_turb_floor ? 1.0e-2 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "K_TURB_CUTOFF_REF",
                MeltingConfig(;
                    [
                        f => (f === :k_turb_cutoff ? 2.0e6 : getfield(m0, f)) for
                        f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "T_SOL",
                MeltingConfig(;
                    [
                        f => (
                            if f === :T_solidus
                                SVector{3,Float64}(1400.0, 1410.0, NaN)
                            else
                                getfield(m0, f)
                            end
                        ) for f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
            (
                "T_LIQ",
                MeltingConfig(;
                    [
                        f => (
                            if f === :T_liquidus
                                SVector{3,Float64}(1800.0, 1810.0, NaN)
                            else
                                getfield(m0, f)
                            end
                        ) for f in fieldnames(MeltingConfig)
                    ]...,
                ),
            ),
        )
        for (needle, m_bad) in melt_mutants
            cfg_b = with_sec(:melting, m_bad)
            ok_b, errs_b = MOSV.check_reference_literals_match_config(cfg_b)
            @test !ok_b
            @test any(occursin(needle, e) for e in errs_b)
            @test MOSV.canonical_config_sha256(cfg_b) != sha
        end

        s0 = base.solver
        for (needle, s_bad) in (
            (
                "ETA_MIN_REF",
                SolverConfig(;
                    seed=s0.seed, hydromech_solver=s0.hydromech_solver, etamin=1.0e13
                ),
            ),
            (
                "ETA_MAX_REF",
                SolverConfig(;
                    seed=s0.seed, hydromech_solver=s0.hydromech_solver, etamax=1.0e24
                ),
            ),
        )
            cfg_b = with_sec(:solver, s_bad)
            ok_b, errs_b = MOSV.check_reference_literals_match_config(cfg_b)
            @test !ok_b
            @test any(occursin(needle, e) for e in errs_b)
            @test MOSV.canonical_config_sha256(cfg_b) != sha
        end

        mat0 = base.materials
        cfg_mm = with_sec(
            :materials,
            MaterialConfig(;
                [
                    f => (
                        if f === :etasolidmm
                            SVector{3,Float64}(2.0e19, 1.0e19, 1.0e18)
                        else
                            getfield(mat0, f)
                        end
                    ) for f in fieldnames(MaterialConfig)
                ]...,
            ),
        )
        ok_mm, errs_mm = MOSV.check_reference_literals_match_config(cfg_mm)
        @test !ok_mm
        @test any(occursin("ETA_SOLID_REF", e) for e in errs_mm)
        @test MOSV.canonical_config_sha256(cfg_mm) != sha

        th0 = base.thermodynamics
        cfg_th = with_sec(
            :thermodynamics,
            ThermalConfig(;
                [
                    f => (f === :tmsolidphase ? 1450.0 : getfield(th0, f)) for
                    f in fieldnames(ThermalConfig)
                ]...,
            ),
        )
        ok_th, errs_th = MOSV.check_reference_literals_match_config(cfg_th)
        @test !ok_th
        @test any(occursin("tmsolidphase", e) for e in errs_th)
        @test MOSV.canonical_config_sha256(cfg_th) != sha
    end

    @testset "Time-step lag, face averaging and unbuffered contrast in 1-D" begin
        lag = MOSV.lag_shift()
        @test isapprox(lag.dT_centre, 1.24; atol=0.02)
        @test 0.0 < lag.dT_centre <= MOSV.metric_limit("centre_tmax_nx65").atol / 3
        @test lag.dT_centre <= MOSV.metric_limit("centre_tmax_nx129").atol / 3
        @test isapprox(lag.dr_mush, 20.0; atol=1e-9)
        @test isapprox(lag.dr_melt, 20.0; atol=1e-9)
        e_c, e_f = lag.coarse.enthalpy_residual, lag.fine.enthalpy_residual
        @test isapprox(e_c, 0.00251; atol=5e-5)
        @test 0.05 < e_f / e_c < 0.15
        @test e_c <= MOSV.ENTHALPY_RESIDUAL_LIMIT
        @test_throws DomainError MOSV.lag_shift(; divisor=0, coarse=lag.coarse)

        # Arithmetic vs harmonic face conductivity in 1-D (Nr = 2000 and coarse Nr = 54)
        fshift = MOSV.face_averaging_shift(; reference=lag.coarse)
        @test isapprox(fshift.dT_centre, -0.0105; atol=0.002)
        @test abs(fshift.dr_mush) <= 10.0
        @test abs(fshift.dr_melt) <= 10.0
        f54 = MOSV.face_averaging_shift(; Nr=54)
        @test isapprox(f54.dT_centre, -2.95; atol=0.02)

        # Thermostat ramp-window vs exponent sensitivity in 1-D
        solve_1d_var(; kw_st...) = MOSV.PSV.solve_porous_sphere(;
            R=MOSV.PSV.R_PLANET,
            phi_init=MOSV.PHI_INIT_CASE,
            t_form_yr=MOSV.T_FORM_YR_CASE,
            t_end_yr=MOSV.N_STEPS_CASE * MOSV.DT_YR_CASE,
            dt_yr=MOSV.DT_YR_CASE,
            Nr=MOSV.NR_REFERENCE,
            stop_cold=false,
            melt=MOSV.MSV.reference_melt(),
            soft_turb=MOSV.reference_soft_turb(; kw_st...),
        )
        tc_nom = MOSV.reference_values(lag.coarse, 65).centre_tmax
        tc_b27 =
            MOSV.reference_values(solve_1d_var(; turb_exponent=2.0 / 7.0), 65).centre_tmax
        tc_f24 =
            MOSV.reference_values(solve_1d_var(; F_start=0.20, F_end=0.40), 65).centre_tmax
        tc_f46 =
            MOSV.reference_values(solve_1d_var(; F_start=0.40, F_end=0.60), 65).centre_tmax
        @test isapprox(tc_b27 - tc_nom, 4.17; atol=0.02)
        @test isapprox(tc_f24 - tc_nom, -41.34; atol=0.02)
        @test isapprox(tc_f46 - tc_nom, 41.46; atol=0.02)

        # Running cell-wise maximum T_max finishes peaking before 1.40 Myr
        ref_25 = MOSV.PSV.solve_porous_sphere(;
            R=MOSV.PSV.R_PLANET,
            phi_init=MOSV.PHI_INIT_CASE,
            t_form_yr=MOSV.T_FORM_YR_CASE,
            t_end_yr=2.50e6,
            dt_yr=MOSV.DT_YR_CASE,
            Nr=MOSV.NR_REFERENCE,
            stop_cold=false,
            melt=MOSV.MSV.reference_melt(),
            soft_turb=MOSV.reference_soft_turb(),
        )
        @test isapprox(lag.coarse.T_max, ref_25.T_max; rtol=1.0e-15, atol=0.0)
        @test isapprox(lag.coarse.t_peak_yr / 1.0e6, 0.91; atol=0.01)
        @test isapprox(lag.coarse.k_max, 705.0; atol=0.5)

        unbuf_14 = MOSV.PSV.solve_porous_sphere(;
            R=MOSV.PSV.R_PLANET,
            phi_init=MOSV.PHI_INIT_CASE,
            t_form_yr=MOSV.T_FORM_YR_CASE,
            t_end_yr=1.40e6,
            dt_yr=MOSV.DT_YR_CASE,
            Nr=MOSV.NR_REFERENCE,
            stop_cold=false,
            melt=MOSV.MSV.reference_melt(),
            soft_turb=MOSV.reference_soft_turb(; enabled=false),
        )
        unbuf_35 = MOSV.PSV.solve_porous_sphere(;
            R=MOSV.PSV.R_PLANET,
            phi_init=MOSV.PHI_INIT_CASE,
            t_form_yr=MOSV.T_FORM_YR_CASE,
            t_end_yr=3.50e6,
            dt_yr=MOSV.DT_YR_CASE,
            Nr=MOSV.NR_REFERENCE,
            stop_cold=false,
            melt=MOSV.MSV.reference_melt(),
            soft_turb=MOSV.reference_soft_turb(; enabled=false),
        )
        v_buf = MOSV.reference_values(lag.coarse, 65)
        v_off14 = MOSV.reference_values(unbuf_14, 65)
        v_off35 = MOSV.reference_values(unbuf_35, 65)
        @test isapprox(v_off14.centre_tmax, 1844.6; atol=0.1)
        @test isapprox(v_off14.centre_tmax - v_buf.centre_tmax, 264.0; atol=0.1)
        @test isapprox(v_off14.plateau_drop, 106.8; atol=0.2)
        @test isapprox(v_off35.centre_tmax, 1984.4; atol=0.2)
        @test isapprox(v_off35.centre_tmax - v_buf.centre_tmax, 403.8; atol=0.2)
        @test isapprox(unbuf_35.t_peak_yr / 1.0e6, 1.98; atol=0.02)
    end

    @testset "1-D planetesimal regime comparison across radius and formation time" begin
        # Nr = 200 with stop_cold = true checks the 4 radii at early and late formation
        reg(R_km, t_Myr) = MOSV.PSV.solve_porous_sphere(;
            R=R_km * 1000.0,
            phi_init=0.25,
            t_form_yr=t_Myr * 1.0e6,
            t_end_yr=3.5e6,
            dt_yr=1.0e4,
            Nr=200,
            stop_cold=true,
            melt=MOSV.MSV.reference_melt(),
            soft_turb=MOSV.reference_soft_turb(),
        )
        r10_05 = reg(10.0, 0.5)
        r10_10 = reg(10.0, 1.0)
        r10_15 = reg(10.0, 1.5)
        r10_25 = reg(10.0, 2.5)
        r20_05 = reg(20.0, 0.5)
        r20_10 = reg(20.0, 1.0)
        r20_15 = reg(20.0, 1.5)
        r20_25 = reg(20.0, 2.5)
        r50_05 = reg(50.0, 0.5)
        r50_10 = reg(50.0, 1.0)
        r50_15 = reg(50.0, 1.5)
        r50_25 = reg(50.0, 2.5)
        r100_05 = reg(100.0, 0.5)
        r100_10 = reg(100.0, 1.0)
        r100_15 = reg(100.0, 1.5)
        r100_25 = reg(100.0, 2.5)

        @test isapprox(r10_05.T_max[1], 1570.9; atol=0.2)
        @test isapprox(r10_10.T_max[1], 1398.8; atol=0.2)
        @test isapprox(r10_15.T_max[1], 1097.5; atol=0.2)
        @test isapprox(r10_25.T_max[1], 757.6; atol=0.2)
        @test iszero(r10_10.F_max[1])
        @test isapprox(r20_05.T_max[1], 1603.6; atol=0.2)
        @test isapprox(r20_10.T_max[1], 1575.5; atol=0.2)
        @test isapprox(r20_15.T_max[1], 1399.7; atol=0.2)
        @test isapprox(r20_25.T_max[1], 801.3; atol=0.2)
        @test isapprox(r50_05.T_max[1], 1679.1; atol=0.2)
        @test isapprox(r50_10.T_max[1], 1599.6; atol=0.2)
        @test isapprox(r50_15.T_max[1], 1530.2; atol=0.2)
        @test isapprox(r50_25.T_max[1], 810.2; atol=0.2)
        @test isapprox(r100_05.T_max[1], 1721.7; atol=0.2)
        @test isapprox(r100_10.T_max[1], 1647.8; atol=0.2)
        @test isapprox(r100_15.T_max[1], 1530.3; atol=0.2)
        @test isapprox(r100_25.T_max[1], 810.2; atol=0.2)
    end

    @testset "Mush radius, plateau drop and metric limit contracts" begin
        R_p, dr = 2.0e4, 500.0
        pts = [
            (r * cos(a), r * sin(a)) for r in (dr / 8):(dr / 4):R_p for
            a in range(0, 2pi; length=65)[1:64]
        ]
        xm, ym = first.(pts), last.(pts)
        rr = hypot.(xm, ym)
        # Step profile crossing T_MUSH = 1520 K at 12.5 km and T_SOL = 1400 K at 15.0 km
        T_step = [r < 12_500.0 ? 1580.0 : (r < 15_000.0 ? 1460.0 : 900.0) for r in rr]
        @test isapprox(
            MOSV.mush_radius(xm, ym, T_step; xc=0.0, yc=0.0, R_p=R_p, dr=dr),
            12_500.0;
            atol=1e-9,
        )
        rc_1d = [(i - 0.5) * R_p / 2000 for i in 1:2000]
        T_1d = [r < 15640.0 ? 1580.0 : (r < 16440.0 ? 1460.0 : 900.0) for r in rc_1d]
        @test isapprox(MOSV.mush_radius_1d(T_1d, R_p), 15640.0; atol=1e-9)

        # Plateau drop: flat vs quadratic vs inverted
        rc_bins = collect(250.0:500.0:19_750.0)
        @test iszero(
            MOSV.plateau_drop(rc_bins, fill(1580.0, length(rc_bins)); r_outer=1.0e4)
        )
        quad = [1600.0 - 100.0 * (r / 1.0e4)^2 for r in rc_bins]
        idx_q = argmin(abs.(rc_bins .- 1.0e4))
        @test isapprox(
            MOSV.plateau_drop(rc_bins, quad; r_outer=1.0e4),
            quad[1] - quad[idx_q];
            rtol=1e-13,
        )
        @test 90.0 < MOSV.plateau_drop(rc_bins, quad; r_outer=1.0e4) < 105.0
        inv_prof = [1500.0 + 20.0 * (r / 1.0e4) for r in rc_bins]
        @test -25.0 < MOSV.plateau_drop(rc_bins, inv_prof; r_outer=1.0e4) < -10.0
        @test !MOSV.metric_passes("plateau_drop_nx65", -0.5, 0.0)
        @test MOSV.metric_passes("plateau_drop_nx65", 4.74, 0.0)
        @test !MOSV.metric_passes("plateau_drop_nx65", 8.51, 0.0)

        @test_throws ArgumentError MOSV.plateau_drop([100.0], [1500.0]; r_outer=100.0)
        @test_throws ArgumentError MOSV.plateau_drop(
            [100.0, 200.0], [1500.0]; r_outer=150.0
        )
        @test_throws ArgumentError MOSV.plateau_drop(
            [300.0, 200.0], [1500.0, 1490.0]; r_outer=250.0
        )
        @test_throws DomainError MOSV.plateau_drop(
            rc_bins, fill(1580.0, length(rc_bins)); r_outer=0.0
        )
        @test_throws DomainError MOSV.plateau_drop(
            rc_bins, fill(1580.0, length(rc_bins)); r_outer=25_000.0
        )

        # Metric limits, units and grid parsing
        @test MOSV.metric_nx("mush_radius_nx65") == 65
        @test MOSV.metric_nx("plateau_drop_nx129") == 129
        @test iszero(MOSV.metric_nx("enthalpy_residual_1d"))
        @test MOSV.metric_units("centre_tmax_nx65") == "K"
        @test MOSV.metric_units("plateau_drop_nx65") == "K"
        @test MOSV.metric_units("mush_radius_nx65") == "m"
        @test MOSV.metric_units("melt_radius_nx65") == "m"
        @test MOSV.metric_units("tmax_rel_l2_hot_nx65") == "dimensionless"
        @test_throws ArgumentError MOSV.metric_units("other_nx65")
        @test isapprox(MOSV.metric_limit("centre_tmax_nx65").atol, 8.0; rtol=1e-14)
        @test isapprox(MOSV.metric_limit("centre_tmax_nx129").atol, 5.0; rtol=1e-14)
        @test isapprox(MOSV.metric_limit("plateau_drop_nx65").atol, 8.0; rtol=1e-14)
        @test isapprox(MOSV.metric_limit("plateau_drop_nx129").atol, 6.0; rtol=1e-14)
        @test MOSV.metric_limit("plateau_drop_nx65").comparison == "le_atol"
        @test isapprox(MOSV.metric_limit("mush_radius_nx65").atol, 1.5 * 750.0; rtol=1e-14)
        @test isapprox(MOSV.metric_limit("melt_radius_nx129").atol, 1.5 * 375.0; rtol=1e-14)
        @test_throws ArgumentError MOSV.metric_limit("centre_tmax_nx33")
    end

    @testset "Soft turbulence off and wrong exponents fail at Nx = 41" begin
        ref = MOSV.reference_solution()
        kw = (;
            Nx=41,
            phi_init=MOSV.PHI_INIT_CASE,
            t_form_yr=MOSV.T_FORM_YR_CASE,
            n_steps=MOSV.N_STEPS_CASE,
            dt_yr=MOSV.DT_YR_CASE,
            reference=ref,
        )
        nominal = MOSV.run_magma_ocean_case(; kw...)
        no_turb = MOSV.run_magma_ocean_case(; kw..., soft_turbulence=false)
        low_exp = MOSV.run_magma_ocean_case(; kw..., turb_exponent=MOSV.TURB_EXP_LOW_MUTANT)
        high_exp = MOSV.run_magma_ocean_case(;
            kw..., turb_exponent=MOSV.TURB_EXP_HIGH_MUTANT
        )
        runs = (nominal, no_turb, low_exp, high_exp)

        pass_c(r) = MOSV.metric_passes("centre_tmax_nx65", r.T_centre, r.T_centre_ref)
        pass_p(r) = MOSV.metric_passes("plateau_drop_nx65", r.plateau_drop, 0.0)
        @test pass_c(nominal) && pass_p(nominal)
        @test !pass_c(no_turb) && !pass_p(no_turb)
        @test !pass_c(low_exp) && !pass_p(low_exp)
        @test !pass_c(high_exp)

        @test isapprox(nominal.T_centre_error, -2.82; atol=0.05)
        @test isapprox(no_turb.T_centre_error, 255.64; atol=0.1)
        @test isapprox(low_exp.T_centre_error, 14.22; atol=0.05)
        @test isapprox(high_exp.T_centre_error, -12.07; atol=0.05)
        @test isapprox(nominal.plateau_drop, 5.28; atol=0.05)
        @test isapprox(no_turb.plateau_drop, 121.54; atol=0.1)
        @test isapprox(low_exp.plateau_drop, 9.11; atol=0.05)
        @test isapprox(high_exp.plateau_drop, 3.89; atol=0.05)
        @test isapprox(nominal.tmax_rel_l2, 0.0314; atol=2e-4)
        @test isapprox(no_turb.tmax_rel_l2, 0.1422; atol=2e-4)
        @test all(r.energy_lse_residual <= 1.0e-9 for r in runs)
        @test all(r.T_min_grid >= MOSV.PSV.T_SPACE - 1.0e-9 for r in runs)
    end

    @testset "Suite record at Nx = 129: static check and pinned 1-D values" begin
        ok, errs = MOSV.check_suite_record_static(MO_SUITE_RECORD)
        @test ok
        @test isempty(errs)
        raw = TOML.parsefile(MO_SUITE_RECORD)
        rec = Dict(String(m["id"]) => m for m in raw["metrics"])
        dT = rec["centre_tmax_nx129"]["value"] - rec["centre_tmax_nx129"]["reference"]
        @test isapprox(dT, -1.63; atol=0.02)
        @test isapprox(rec["plateau_drop_nx129"]["value"], 4.29; atol=0.02)
        @test isapprox(rec["mush_radius_nx129"]["value"], 15187.5; atol=1e-9)
        @test isapprox(rec["mush_radius_nx129"]["reference"], 15640.0; atol=1e-9)
        @test isapprox((15187.5 - 15640.0) / 375.0, -1.21; atol=0.01)
        @test isapprox(rec["melt_radius_nx129"]["value"], 16031.25; atol=1e-9)
        @test isapprox(rec["melt_radius_nx129"]["reference"], 16440.0; atol=1e-9)
        @test isapprox((16031.25 - 16440.0) / 375.0, -1.09; atol=0.01)
        @test isapprox(MOSV.metric_limit("mush_radius_nx129").atol, 562.5; rtol=1e-14)
        @test isapprox(rec["tmax_rel_l2_hot_nx129"]["value"], 0.0188; atol=1e-4)

        cli = _mo_capture(() -> run_validation_cli(["--check", MO_SUITE_RECORD]))
        @test cli.value == 0
        @test static_validation_case("magma_ocean_conductivity_suite").record ==
            DEFAULT_MAGMA_OCEAN_CONDUCTIVITY_SUITE_RECORD

        mktempdir() do dir
            for (mid, changes, needle) in (
                ("centre_tmax_nx129", ("atol" => 15.0,), "differs from the limit"),
                (
                    "centre_tmax_nx129",
                    ("reference" => 1550.0, "value" => 1550.0),
                    "pinned 1-D value",
                ),
                ("plateau_drop_nx129", ("value" => 25.0,), "fails comparison"),
                ("mush_radius_nx129", ("units" => "km",), "units"),
            )
                bad = deepcopy(raw)
                i = findfirst(m -> m["id"] == mid, bad["metrics"])
                for (k, v) in changes
                    bad["metrics"][i][k] = v
                end
                ok_b, errs_b = MOSV.check_suite_record_static(
                    _mo_write_toml(joinpath(dir, "bad.toml"), bad); check_doc_sync=false
                )
                @test !ok_b
                @test any(occursin(needle, e) for e in errs_b)
            end
        end
        @test !first(MOSV.check_suite_record_static(MO_SUITE_RECORD * ".missing"))
        @test_throws DomainError MOSV.check_suite_record_static("")
    end

    bench = MOSV.evaluate_magma_ocean_benchmark()

    @testset "CI run at Nx = 65 against the 1-D reference" begin
        r = bench.runs[65]
        @test isapprox(r.T_centre_ref, 1580.64; atol=0.02)
        @test isapprox(r.T_centre, 1578.38; atol=0.02)
        @test isapprox(r.T_centre_error, -2.26; atol=0.02)
        @test isapprox(r.plateau_drop, 4.74; atol=0.02)
        @test isapprox(r.plateau_drop_ref, 3.33; atol=0.02)
        @test isapprox(r.r_mush, 15000.0; atol=1e-9)
        @test isapprox(r.r_mush_ref, 15640.0; atol=1e-9)
        @test isapprox(r.r_mush - r.r_mush_ref, -640.0; atol=1e-9)
        @test isapprox(r.r_mush_error_cells, -0.85; atol=0.01)
        @test isapprox(r.r_melt, 15750.0; atol=1e-9)
        @test isapprox(r.r_melt_ref, 16440.0; atol=1e-9)
        @test isapprox(r.r_melt - r.r_melt_ref, -690.0; atol=1e-9)
        @test isapprox(r.r_melt_error_cells, -0.92; atol=0.01)
        @test isapprox(r.tmax_rel_l2, 0.0283; atol=1e-4)
        @test r.energy_lse_residual <= 1.0e-9
        @test r.T_min_grid >= MOSV.PSV.T_SPACE - 1.0e-9
        @test isapprox(bench.metrics["enthalpy_residual_1d"], 0.00251; atol=5e-5)

        v129 = MOSV.reference_values(bench.reference, 129)
        @test isapprox(v129.plateau_drop, 3.42; atol=0.01)
        @test isapprox(
            v129.centre_tmax, MOSV.suite_reference_value("centre_tmax_nx129"); rtol=1e-12
        )
        @test isapprox(
            v129.mush_radius, MOSV.suite_reference_value("mush_radius_nx129"); rtol=1e-12
        )
        @test isapprox(
            v129.melt_radius, MOSV.suite_reference_value("melt_radius_nx129"); rtol=1e-12
        )
    end

    @testset "CI record: reproducibility, tampering and export round trip" begin
        ok, errs = MOSV.check_validation_reproducibility(bench, MO_RECORD)
        @test ok
        @test isempty(errs)
        @test run_magma_ocean_conductivity_validation_check(MO_RECORD; bench=bench)
        @test run_validation_cli(
            ["--check", MO_RECORD]; magma_ocean_conductivity_bench=bench
        ) == 0
        @test last(VALIDATION_CASES).id == "magma_ocean_conductivity"

        raw = TOML.parsefile(MO_RECORD)
        @test [String(m["id"]) for m in raw["metrics"]] == collect(MOSV.EXPECTED_METRIC_IDS)
        @test [Int(m["nx"]) for m in raw["metrics"]] == [65, 65, 65, 65, 65, 65, 0]

        mktempdir() do dir
            i_c = findfirst(m -> m["id"] == "centre_tmax_nx65", raw["metrics"])
            wide = deepcopy(raw)
            wide["metrics"][i_c]["atol"] = 20.0
            ok_w, errs_w = MOSV.check_validation_reproducibility(
                bench, _mo_write_toml(joinpath(dir, "wide.toml"), wide)
            )
            @test !ok_w
            @test any(occursin("differs from the limit of the code", e) for e in errs_w)

            moved = deepcopy(raw)
            moved["metrics"][i_c]["value"] += 3.0
            moved["metadata"]["julia_version"] = string(VERSION)
            ok_m, errs_m = MOSV.check_validation_reproducibility(
                bench, _mo_write_toml(joinpath(dir, "moved.toml"), moved)
            )
            @test !ok_m
            @test any(occursin("differs from recorded value", e) for e in errs_m)

            # A +5.0 K thermostat shift passes the 8.0 K limit, caught by 0.00906 K band
            @test isapprox(MOSV.SEED_BAND_CENTRE_TMAX, 0.00906; rtol=1e-14)
            @test isapprox(MOSV.SEED_BAND_PLATEAU_DROP, 0.0237; rtol=1e-14)
            hot_metrics = copy(bench.metrics)
            hot_metrics["centre_tmax_nx65"] += 5.0
            hot = (; bench..., metrics=hot_metrics)
            ok_h, errs_h = MOSV.check_validation_reproducibility(hot, MO_RECORD)
            @test !ok_h
            @test !any(occursin("fails its limit", e) for e in errs_h)
            @test any(occursin("physics regression band", e) for e in errs_h)

            low_metrics = copy(bench.metrics)
            low_metrics["centre_tmax_nx65"] -= 12.0
            low = (; bench..., metrics=low_metrics)
            ok_l, errs_l = MOSV.check_validation_reproducibility(low, MO_RECORD)
            @test !ok_l
            @test any(occursin("fails its limit", e) for e in errs_l)
            @test any(occursin("physics regression band", e) for e in errs_l)

            out = joinpath(dir, "exported.toml")
            commit = String(raw["metadata"]["src_commit"])
            export_cli(args, b) = MagmaOceanConductivityExport.run_export_cli(
                args; src_commit=commit, bench=b, check_clean_tree=false
            )
            @test _mo_capture(() -> export_cli([out], bench)).value == 0
            @test first(MOSV.check_validation_reproducibility(bench, out))
            @test _mo_capture(() -> export_cli(["--bad", out], bench)).value == 2
            refused = _mo_capture(() -> export_cli([out], low))
            @test refused.value == 1
            @test occursin("accept-physics-change", refused.err)
            refused_band = _mo_capture(() -> export_cli([out], hot))
            @test refused_band.value == 1
            @test occursin("accept-physics-change", refused_band.err)
            @test _mo_capture(() -> export_cli(["--accept-physics-change", out], hot)).value ==
                0
            @test_throws DomainError MagmaOceanConductivityExport.run_export_cli(
                [out]; src_commit="abc", bench=bench, check_clean_tree=false
            )
            @test_throws ArgumentError MagmaOceanConductivityExport.metric_description(
                "unknown_metric_nx65"
            )
            @test_throws ErrorException MagmaOceanConductivityExport.magma_ocean_validation_module((
                Base,
            ))

            # Suite export round-trip and CSV validation
            sdir = mkpath(joinpath(dir, "suite"))
            write(joinpath(sdir, "header.csv"), join(MOSV.SUITE_RESULT_HEADER, ",") * "\n")
            s_rec = Dict(
                String(m["id"]) => m for m in TOML.parsefile(MO_SUITE_RECORD)["metrics"]
            )
            s_row = (
                "nx129",
                129,
                "true",
                MOSV.TURB_EXP_REF,
                MOSV.PSV.L_BOX / 128,
                s_rec["centre_tmax_nx129"]["value"],
                MOSV.suite_reference_value("centre_tmax_nx129"),
                s_rec["plateau_drop_nx129"]["value"],
                3.42237,
                s_rec["mush_radius_nx129"]["value"],
                MOSV.suite_reference_value("mush_radius_nx129"),
                s_rec["melt_radius_nx129"]["value"],
                MOSV.suite_reference_value("melt_radius_nx129"),
                s_rec["tmax_rel_l2_hot_nx129"]["value"],
                2.0e-13,
                0.002513,
                374.8,
            )
            write(joinpath(sdir, "result_nx129.csv"), join(s_row, ",") * "\n")
            sout = joinpath(dir, "exported_suite.toml")
            @test _mo_capture(() -> export_cli(["--suite", sdir, sout], nothing)).value == 0
            @test first(MOSV.check_suite_record_static(sout; check_doc_sync=false))

            bad_sdir = mkpath(joinpath(dir, "bad_suite"))
            @test _mo_capture(() -> export_cli(["--suite", bad_sdir, sout], nothing)).value ==
                1
            write(joinpath(bad_sdir, "header.csv"), "bad,header\n")
            @test _mo_capture(() -> export_cli(["--suite", bad_sdir, sout], nothing)).value ==
                1

            # Dirty source tree refusal check
            @test isapprox(
                MOSV.metric_limit("tmax_rel_l2_hot_nx65").atol, 0.04; atol=1.0e-12
            )
            @test isapprox(
                MOSV.metric_limit("tmax_rel_l2_hot_nx129").atol, 0.025; atol=1.0e-12
            )
            fake_repo = mkpath(joinpath(dir, "fake_repo"))
            run(`git -C $fake_repo init -q`)
            mkpath(joinpath(fake_repo, "src"))
            write(joinpath(fake_repo, "src", "a.jl"), "x = 1\n")
            run(`git -C $fake_repo add src/a.jl`)
            run(`git -C $fake_repo -c user.name=t -c user.email=t@t commit -q -m init`)
            write(joinpath(fake_repo, "src", "a.jl"), "x = 2\n")
            dirty_run = _mo_capture(
                () -> MagmaOceanConductivityExport.run_export_cli(
                    [joinpath(fake_repo, "exported.toml")];
                    src_commit=commit,
                    bench=bench,
                    check_clean_tree=true,
                    repo_root=fake_repo,
                ),
            )
            @test dirty_run.value == 1
            @test occursin("Refusing to export validation record", dirty_run.err)
        end
        @test !first(MOSV.check_validation_reproducibility(bench, MO_RECORD * ".x"))
        @test_throws DomainError MOSV.check_validation_reproducibility(bench, "")
    end
end
