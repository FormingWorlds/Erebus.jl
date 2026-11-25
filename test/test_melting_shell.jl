using Test
using TOML
using Erebus
isdefined(@__MODULE__, :MeltingShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_melting_shell.jl"))
isdefined(@__MODULE__, :run_validation_cli) ||
    include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
isdefined(@__MODULE__, :MeltingShellExport) ||
    include(joinpath(@__DIR__, "..", "tools", "export_validation_melting_shell.jl"))
const MSV = MeltingShellValidation
const MS_RECORD = joinpath(@__DIR__, "..", "benchmarks", "records", "melting_shell.toml")
const MS_SUITE_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "melting_shell_suite.toml"
)

"""
Write the record dictionary `dict` as TOML to `path` and return `path`.
"""
function _ms_write_toml(path::String, dict::AbstractDict)
    open(io -> TOML.print(io, dict), path, "w")
    return path
end

"""
Call `f()` with standard output and error sent to files; return the value and the text.
"""
function _ms_capture(f)
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
Set every rock marker of `st` to temperature `T` and dry porosity `phi_d`, then run the
heating and the marker-to-grid step of `cfg`. Returns the rock marker indices.
"""
function _ms_marker_step!(st, coords, ws, cfg, T, phi_d)
    mk = st.markers
    rock = findall(<(3), mk.tm)
    mk.tkm[rock] .= T
    mk.phi_d[rock] .= phi_d
    Erebus.radiogenic_heating!(st, coords, cfg)
    Erebus.interpolate_markers_to_grid!(
        st, coords, cfg; thread_buffers=ws.thread_buffers, interp_arrays=ws.interp_arrays
    )
    return rock
end

"""
Initialize a state for `cfg` and close the telemetry stream.
"""
function _ms_state(cfg)
    st, coords, ws, tel, _ = Erebus.init_simulation(cfg; output_path=mktempdir())
    tel !== nothing && close(tel)
    return st, coords, ws
end

@testset "Melting Shell Validation" begin
    @testset "Marker heat capacity with latent heat at unsintered porosity" begin
        rhocp_s, rho_s, dT = 3.3e6, 3300.0, 400.0
        # At L 4e5 the latent term rho L / dT equals rho cp; at 6e5 it is 4.95e6, so this
        # loop tells the 2 terms apart
        for L in (4.0e5, 6.0e5)
            cfg = MSV.make_melting_shell_config(;
                Nx=17, phi_init=0.25, t_form_yr=1.3e6, n_steps=1, dt_yr=5.0e4, L_melt=L
            )
            @test cfg.melting.active
            @test isapprox(cfg.melting.rho_melt, rho_s; rtol=1e-14)
            st, coords, ws = _ms_state(cfg)
            # Inside the interval: (1 - phi_d) (rho cp + rho L / (T_l - T_s))
            rock = _ms_marker_step!(st, coords, ws, cfg, 1600.0, 0.40)
            expected = 0.6 * (rhocp_s + rho_s * L / dT)
            @test all(isapprox.(st.markers.rhocptotalm[rock], expected; rtol=1e-12))
            # At the interval edges the latent term is off (strict interval)
            for T in (1400.0, 1800.0)
                rock = _ms_marker_step!(st, coords, ws, cfg, T, 0.40)
                @test all(
                    isapprox.(st.markers.rhocptotalm[rock], 0.6 * rhocp_s; rtol=1e-12)
                )
            end
        end
    end

    @testset "Melting leaves the radiogenic heat of a marker unchanged" begin
        kw = (; Nx=17, phi_init=0.25, t_form_yr=1.3e6, n_steps=1, dt_yr=5.0e4)
        cfg_on = MSV.make_melting_shell_config(; kw...)
        cfg_off = MSV.PorousShellValidation.make_porous_shell_config(; kw...)
        hr = map((cfg_on, cfg_off)) do cfg
            st, coords, ws = _ms_state(cfg)
            rock = _ms_marker_step!(st, coords, ws, cfg, 1600.0, 0.40)
            st.markers.hrtotalm[rock]
        end
        @test all(>(0.0), hr[1])
        # Zero tolerances: the radiogenic heat agrees exactly
        @test isapprox(hr[1], hr[2]; rtol=0.0, atol=0.0)
    end

    @testset "Reference heat capacity matches the Erebus law" begin
        # solve_porous_sphere takes its heat capacity from this law (PSV.rhocp_app_ref)
        for L in (4.0e5, 6.0e5),
            phi in (0.0, 0.4),
            T in (1000.0, 1400.0, 1401.0, 1799.0, 1800.0, 2000.0)

            src = Erebus.rhocp_apparent_silicate(
                T, 0.0, 3.3e6, 3300.0, 1; T_sol=1400.0, T_liq=1800.0, L_melt=L
            )
            @test isapprox(MSV.rhocp_app_ref(T, phi; L_melt=L), (1 - phi) * src; rtol=1e-14)
        end
        @test iszero(MSV.melt_fraction_ref(1300.0))
        @test isapprox(MSV.melt_fraction_ref(1500.0), 0.25; atol=1e-14)
        @test isapprox(MSV.melt_fraction_ref(1900.0), 1.0; atol=1e-14)
    end

    @testset "1-D reference: L = 0 equals the run without melting" begin
        kw = (;
            R=2.0e4, phi_init=0.25, t_form_yr=1.3e6, t_end_yr=2.0e6, dt_yr=5.0e4, Nr=200
        )
        a = MSV.PSV.solve_porous_sphere(; kw...)
        b = MSV.PSV.solve_porous_sphere(;
            kw..., melt=(; L_melt=0.0, T_sol=1400.0, T_liq=1800.0)
        )
        # Zero tolerances: the 2 runs agree exactly
        @test isapprox(a.T, b.T; rtol=0.0, atol=0.0)
        @test isapprox(a.T_max, b.T_max; rtol=0.0, atol=0.0)
        @test isapprox(a.phi, b.phi; rtol=0.0, atol=0.0)
        @test isapprox(a.d_shell, b.d_shell; rtol=0.0, atol=0.0)
        # Without latent heat the stepwise enthalpy closes to round-off
        @test a.enthalpy_residual < 1e-10
    end

    @testset "1-D reference: adiabatic centre brackets the enthalpy of the deposited heat" begin
        # R 200 km: conduction reaches ~5 km in 1 Myr, so the centre cell is adiabatic
        c_s, c_lat = 3.3e6, 3300.0 * 4.0e5 / 400.0
        t_form, t_end = 1.1e6, 1.0e6
        errs = Float64[]
        for k in 0:4
            dt = 5.0e4 / 2^k
            res = MSV.PSV.solve_porous_sphere(;
                R=2.0e5,
                phi_init=0.0,
                t_form_yr=t_form,
                t_end_yr=t_end,
                dt_yr=dt,
                Nr=400,
                melt=(; L_melt=4.0e5, T_sol=1400.0, T_liq=1800.0),
                stop_cold=false,
            )
            Tc = res.T[1]
            @test 1400.0 < Tc < 1800.0
            q = 0.0
            for n in 0:(res.n_steps - 1)
                t = (t_form + n * dt) * MSV.PSV.YEARLENGTH_S
                q +=
                    MSV.PSV.RHO_BULK *
                    MSV.PSV.q_al_step(t, dt * MSV.PSV.YEARLENGTH_S, MSV.PSV.RATIO_AL_L16) *
                    dt *
                    MSV.PSV.YEARLENGTH_S
            end
            err = MSV.enthalpy_s_ref(Tc) - q
            # The entry step takes no latent heat: excess between 0 and c_lat dT_step
            dT_step =
                MSV.PSV.RHO_BULK *
                MSV.PSV.q_al_step(
                    t_form * MSV.PSV.YEARLENGTH_S,
                    dt * MSV.PSV.YEARLENGTH_S,
                    MSV.PSV.RATIO_AL_L16,
                ) *
                dt *
                MSV.PSV.YEARLENGTH_S / c_s
            @test -1e-9 * q <= err <= c_lat * dT_step * (1 + 1e-9)
            push!(errs, err)
        end
        @test all(diff(errs) .<= 1e-9 * abs(errs[1]))
        @test errs[end] < errs[1]
    end

    @testset "1-D reference: enthalpy closure with latent heat converges in dt" begin
        run(dt) = MSV.PSV.solve_porous_sphere(;
            R=2.0e4,
            phi_init=0.25,
            t_form_yr=1.3e6,
            t_end_yr=3.0e6,
            dt_yr=dt,
            Nr=200,
            melt=(; L_melt=4.0e5, T_sol=1400.0, T_liq=1800.0),
        )
        r = [run(5.0e4 / 2^k) for k in 0:3]
        @test MSV.T_SOL + 50.0 < r[1].T_peak < MSV.T_LIQ
        e = [x.enthalpy_residual for x in r]
        @test 1e-6 < e[1] < 0.1
        @test all(0.3 .< e[2:end] ./ e[1:(end - 1)] .< 0.7)
        # Stefan number above 1 lowers the peak temperature
        r6 = MSV.PSV.solve_porous_sphere(;
            R=2.0e4,
            phi_init=0.25,
            t_form_yr=1.3e6,
            t_end_yr=3.0e6,
            dt_yr=5.0e4,
            Nr=200,
            melt=(; L_melt=6.0e5, T_sol=1400.0, T_liq=1800.0),
        )
        @test r6.T_peak < r[1].T_peak - 1.0
        @test isapprox(r[1].F_max[1], MSV.melt_fraction_ref(r[1].T_max[1]); atol=1e-14)
    end

    @testset "Melt radius from centre-anchored bins of the peak temperature" begin
        # Polar lattice of markers: 16 sectors, radial spacing dr / 4
        R_p, dr = 2.0e4, 500.0
        lattice(Rmax) = [
            (r * cos(a), r * sin(a)) for r in (dr / 8):(dr / 4):Rmax for
            a in range(0, 2pi; length=97)[1:96]
        ]
        pts = lattice(R_p)
        xm, ym = first.(pts), last.(pts)
        rr = hypot.(xm, ym)
        f(prof) = MSV.melt_radius(xm, ym, prof.(rr); xc=0.0, yc=0.0, R_p=R_p, dr=dr)
        @test isapprox(f(r -> r < 7500.0 ? 1600.0 : 1000.0), 7500.0; atol=1e-9)
        @test iszero(f(r -> 1000.0))
        @test isapprox(f(r -> 1600.0), R_p; atol=1e-9)
        # A hot ring outside a cold bin does not extend the radius
        ring = f(r -> (r < 5000.0 || 6000.0 <= r < 7000.0) ? 1600.0 : 1000.0)
        @test isapprox(ring, 5000.0; atol=1e-9)
        # Exactly at the solidus a bin is not molten
        @test isapprox(f(r -> r < 7500.0 ? 1600.0 : 1400.0), 7500.0; atol=1e-9)
        @test iszero(f(r -> r < 7500.0 ? 1400.0 : 1000.0))
        # Empty bins are skipped
        keep = .!(3000.0 .<= rr .< 3500.0)
        r_keep = MSV.melt_radius(
            xm[keep],
            ym[keep],
            [r < 7500.0 ? 1600.0 : 1000.0 for r in rr[keep]];
            xc=0.0,
            yc=0.0,
            R_p=R_p,
            dr=dr,
        )
        @test isapprox(r_keep, 7500.0; atol=1e-9)
        # Median over sectors: 8 sectors to 7.5 km, 8 sectors to 8.5 km
        ang = mod.(atan.(ym, xm), 2pi)
        T2 = [
            (a < pi ? r < 7500.0 : r < 8500.0) ? 1600.0 : 1000.0 for (r, a) in zip(rr, ang)
        ]
        @test isapprox(
            MSV.melt_radius(xm, ym, T2; xc=0.0, yc=0.0, R_p=R_p, dr=dr), 8000.0; atol=1e-9
        )
        # R_p / dr not an integer: the outer bin ends at R_p; markers at r >= R_p are ignored
        R2 = 2.02e4
        pts2 = lattice(R2 + 400.0)
        x2, y2 = first.(pts2), last.(pts2)
        r2 = hypot.(x2, y2)
        T3 = [r < R2 ? 1600.0 : 1000.0 for r in r2]
        @test isapprox(
            MSV.melt_radius(x2, y2, T3; xc=0.0, yc=0.0, R_p=R2, dr=dr), R2; atol=1e-9
        )
        # Offset centre
        r_off = MSV.melt_radius(
            xm .+ 1.0e4,
            ym .- 3.0e3,
            [r < 7500.0 ? 1600.0 : 1000.0 for r in rr];
            xc=1.0e4,
            yc=-3.0e3,
            R_p=R_p,
            dr=dr,
        )
        @test isapprox(r_off, 7500.0; atol=1e-9)
        @test_throws ArgumentError MSV.melt_radius(
            [0.0], [0.0, 1.0], [1.0]; xc=0.0, yc=0.0, R_p=R_p, dr=dr
        )
        @test_throws DomainError MSV.melt_radius(
            [0.0], [0.0], [1.0]; xc=0.0, yc=0.0, R_p=R_p, dr=0.0
        )
    end

    @testset "Melt radius of a 1-D profile" begin
        R, nr = 2.0e4, 2000
        rc = [(i - 0.5) * R / nr for i in 1:nr]
        @test isapprox(
            MSV.melt_radius_1d([r < 13130.0 ? 1600.0 : 1000.0 for r in rc], R),
            13130.0;
            atol=1e-9,
        )
        # A cell at the solidus is not molten; a hot ring outside does not count
        r_sol = MSV.melt_radius_1d([r < 5000.0 ? 1600.0 : 1400.0 for r in rc], R)
        @test isapprox(r_sol, 5000.0; atol=1e-9)
        ring = [(r < 5000.0 || 8000.0 <= r < 9000.0) ? 1600.0 : 1000.0 for r in rc]
        @test isapprox(MSV.melt_radius_1d(ring, R), 5000.0; atol=1e-9)
        @test iszero(MSV.melt_radius_1d(fill(1000.0, nr), R))
        @test isapprox(MSV.melt_radius_1d(fill(1600.0, nr), R), R; atol=1e-9)
        @test_throws ArgumentError MSV.melt_radius_1d(Float64[], R)
        @test_throws DomainError MSV.melt_radius_1d([1600.0], -1.0)
    end

    @testset "Metric ids, units and limits of the melting shell records" begin
        @test MSV.metric_nx("centre_tmax_nx65") == 65
        @test MSV.metric_nx("melt_radius_nx129") == 129
        @test iszero(MSV.metric_nx("enthalpy_residual_1d"))
        @test_throws ArgumentError MSV.metric_nx("centre_tmax")
        @test MSV.metric_units("centre_tmax_nx65") == "K"
        @test MSV.metric_units("melt_radius_nx65") == "m"
        @test MSV.metric_units("tmax_rel_l2_hot_nx65") == "dimensionless"
        @test MSV.metric_units("enthalpy_residual_1d") == "dimensionless"
        @test_throws ArgumentError MSV.metric_units("shell_depth_nx65")
        # Limits fixed before the records: 20 K and 10 K, 1.5 cells, 0.04 and 0.025
        lim(id) = MSV.metric_limit(id).atol
        @test isapprox(lim("centre_tmax_nx65"), 20.0; rtol=1e-14)
        @test isapprox(lim("centre_tmax_nx129"), 10.0; rtol=1e-14)
        @test isapprox(lim("melt_radius_nx65"), 1.5 * 750.0; rtol=1e-14)
        @test isapprox(lim("melt_radius_nx129"), 1.5 * 375.0; rtol=1e-14)
        @test isapprox(lim("tmax_rel_l2_hot_nx65"), 0.04; rtol=1e-14)
        @test isapprox(lim("tmax_rel_l2_hot_nx129"), 0.025; rtol=1e-14)
        @test isapprox(lim("energy_lse_residual_nx65"), 1.0e-9; rtol=1e-14)
        @test isapprox(lim("enthalpy_residual_1d"), 0.01; rtol=1e-14)
        @test MSV.metric_limit("tmax_rel_l2_hot_nx65").comparison == "le_atol"
        @test MSV.metric_limit("centre_tmax_nx65").comparison == "isapprox"
        @test isapprox(MSV.TMAX_L2_T_MIN, 1100.0; rtol=1e-14)
        @test_throws ArgumentError MSV.metric_limit("centre_tmax_nx33")
        @test_throws ArgumentError MSV.metric_limit("unknown_nx65")
        @test MSV.metric_passes("centre_tmax_nx65", 1645.0, 1661.0)
        @test !MSV.metric_passes("centre_tmax_nx65", 1640.0, 1661.0)
        @test !MSV.metric_passes("tmax_rel_l2_hot_nx65", 0.02, 1.0)
        @test !MSV.metric_passes("tmax_rel_l2_hot_nx65", NaN, 0.0)
    end

    @testset "Reference literals of the melting shell match the configuration" begin
        ok, errs = MSV.check_reference_literals_match_config()
        @test ok
        @test isempty(errs)
        cfg_L = MSV.make_case_config(; L_melt=5.0e5)
        ok_L, errs_L = MSV.check_reference_literals_match_config(cfg_L)
        @test !ok_L
        @test any(occursin("L_MELT", e) for e in errs_L)
        # The hash covers the configuration and the constants of the case
        sha = MSV.canonical_config_sha256(MSV.make_case_config())
        @test length(sha) == 64
        @test MSV.canonical_config_sha256(cfg_L) != sha
        @test MSV.canonical_config_sha256(MSV.make_case_config(; Nx=33)) != sha
        @test MSV.suite_config_sha256() != sha
        text = MSV.canonical_config_text(MSV.make_case_config())
        @test occursin("case.CENTRE_TMAX_ATOL=", text)
        @test occursin("melting.L_melt=400000.0", text)
    end

    @testset "Time lag of the apparent heat capacity in the 1-D reference" begin
        # The reference at 25 kyr against 2.5 kyr: the lag warms the centre by 1.42 K
        # and moves the melt radius out by 50 m, within a third of the 20 K limit.
        lag = MSV.lag_shift()
        @test isapprox(lag.dT_centre, 1.42; atol=0.01)
        @test 0.0 < lag.dT_centre <= MSV.metric_limit("centre_tmax_nx65").atol / 3
        @test lag.dT_centre <= MSV.metric_limit("centre_tmax_nx129").atol / 3
        @test isapprox(lag.dr_melt, 50.0; atol=1e-9)
        @test abs(lag.dr_melt) <= 750.0 / 3
        # The enthalpy residual is first order in the step: 0.0061 at 25 kyr
        e_c, e_f = lag.coarse.enthalpy_residual, lag.fine.enthalpy_residual
        @test isapprox(e_c, 0.0061; atol=5e-5)
        @test 0.05 < e_f / e_c < 0.15
        @test e_c <= MSV.ENTHALPY_RESIDUAL_LIMIT
        @test_throws DomainError MSV.lag_shift(; divisor=0, coarse=lag.coarse)
    end

    @testset "L = 0 and 1.5 L fail, 0.75 L passes the centre limit at Nx = 41" begin
        # Nx 41 is the coarsest grid where the run at L passes the 20 K limit of Nx 65.
        ref = MSV.reference_solution()
        kw = (;
            Nx=41,
            phi_init=MSV.PHI_INIT_CASE,
            t_form_yr=MSV.T_FORM_YR_CASE,
            n_steps=MSV.N_STEPS_CASE,
            dt_yr=MSV.DT_YR_CASE,
            reference=ref,
        )
        L_high = MSV.L_MUTANT_FACTOR * MSV.L_MELT
        nominal = MSV.run_melting_shell_case(; kw..., L_melt=MSV.L_MELT)
        no_lat = MSV.run_melting_shell_case(; kw..., L_melt=0.0)
        high_lat = MSV.run_melting_shell_case(; kw..., L_melt=L_high)
        low_lat = MSV.run_melting_shell_case(; kw..., L_melt=0.75 * MSV.L_MELT)
        runs = (nominal, no_lat, high_lat, low_lat)
        pass(r) = MSV.metric_passes("centre_tmax_nx65", r.T_centre, r.T_centre_ref)
        @test pass(nominal)
        @test !pass(no_lat)
        @test !pass(high_lat)
        # A latent heat 25 % low passes: the grid error sits on the low side of the window
        @test pass(low_lat)
        # Errors against the 1-D model at L: -16.0 K, +89.5 K, -46.2 K and +3.0 K
        @test isapprox(nominal.T_centre_error, -16.0; atol=0.05)
        @test isapprox(no_lat.T_centre_error, 89.5; atol=0.05)
        @test isapprox(high_lat.T_centre_error, -46.2; atol=0.05)
        @test isapprox(low_lat.T_centre_error, 3.0; atol=0.05)
        @test isapprox(20.0 - abs(nominal.T_centre_error), 4.0; atol=0.05)
        # The shift from the latent heat alone, on the same grid, exceeds the limit
        limit = MSV.metric_limit("centre_tmax_nx65").atol
        @test no_lat.T_centre - nominal.T_centre > limit
        @test nominal.T_centre - high_lat.T_centre > limit
        # The shift of 0.75 L is 10 times the seed band of the Nx 65 record
        @test low_lat.T_centre - nominal.T_centre > 10 * MSV.SEED_BAND_CENTRE_TMAX
        # The melt radius stays within its 1.5-cell band for all 4 runs
        @test all(abs(r.r_melt_error_cells) <= 1.5 for r in runs)
        # L2 of the hot bins: 0.0289 at L, 0.0481 at 0, 0.0391 at 1.5 L, 0.0254 at 0.75 L
        @test isapprox(nominal.tmax_rel_l2, 0.0289; atol=1e-4)
        @test isapprox(no_lat.tmax_rel_l2, 0.0481; atol=1e-4)
        @test isapprox(high_lat.tmax_rel_l2, 0.0391; atol=1e-4)
        @test isapprox(low_lat.tmax_rel_l2, 0.0254; atol=1e-4)
        @test all(r.energy_lse_residual <= 1.0e-9 for r in runs)
    end

    @testset "Suite record at Nx = 129: static check and pinned 1-D values" begin
        ok, errs = MSV.check_suite_record_static(MS_SUITE_RECORD)
        @test ok
        @test isempty(errs)
        raw = TOML.parsefile(MS_SUITE_RECORD)
        rec = Dict(String(m["id"]) => m for m in raw["metrics"])
        # Erebus at Nx 129: centre 5.9 K below the 1-D value, melt radius 12750 m
        # against 13130 m (380 m inside), L2 0.0136
        dT = rec["centre_tmax_nx129"]["value"] - rec["centre_tmax_nx129"]["reference"]
        @test isapprox(dT, -5.9; atol=0.05)
        r129, r129_ref = rec["melt_radius_nx129"]["value"],
        rec["melt_radius_nx129"]["reference"]
        @test isapprox(r129, 12750.0; atol=1e-9)
        @test isapprox(r129_ref, 13130.0; atol=1e-9)
        @test isapprox(r129 - r129_ref, -380.0; atol=1e-9)
        @test isapprox(rec["tmax_rel_l2_hot_nx129"]["value"], 0.0136; atol=1e-4)
        cli = _ms_capture(() -> run_validation_cli(["--check", MS_SUITE_RECORD]))
        @test cli.value == 0
        @test static_validation_case("melting_shell_suite").record ==
            DEFAULT_MELTING_SHELL_SUITE_RECORD
        mktempdir() do dir
            # A value outside its recorded tolerance fails the schema first
            for (mid, changes, needle) in (
                ("centre_tmax_nx129", ("atol" => 30.0,), "differs from the limit"),
                (
                    "centre_tmax_nx129",
                    ("reference" => 1600.0, "value" => 1600.0),
                    "pinned 1-D value",
                ),
                ("melt_radius_nx129", ("value" => 10000.0,), "fails comparison"),
                ("tmax_rel_l2_hot_nx129", ("units" => "K",), "units"),
            )
                bad = deepcopy(raw)
                i = findfirst(m -> m["id"] == mid, bad["metrics"])
                for (key, val) in changes
                    bad["metrics"][i][key] = val
                end
                ok_b, errs_b = MSV.check_suite_record_static(
                    _ms_write_toml(joinpath(dir, "bad.toml"), bad); check_doc_sync=false
                )
                @test !ok_b
                @test any(occursin(needle, e) for e in errs_b)
            end
            cfg_bad = deepcopy(raw)
            cfg_bad["config"]["Nx"] = 65
            ok_c, errs_c = MSV.check_suite_record_static(
                _ms_write_toml(joinpath(dir, "cfg.toml"), cfg_bad); check_doc_sync=false
            )
            @test !ok_c
            @test any(occursin("config.Nx", e) for e in errs_c)
        end
        @test !first(MSV.check_suite_record_static(MS_SUITE_RECORD * ".missing"))
        @test_throws DomainError MSV.check_suite_record_static("")
    end

    bench = MSV.evaluate_melting_shell_benchmark()

    @testset "CI run at Nx = 65 against the 1-D reference" begin
        r = bench.runs[65]
        # 1-D centre 1661.3 K; Erebus 11.6 K lower, melt radius 755 m (1.01 cells) inside
        @test isapprox(r.T_centre_ref, 1661.3; atol=0.05)
        @test isapprox(r.T_centre_error, -11.6; atol=0.05)
        @test MSV.T_SOL < r.T_centre < MSV.T_LIQ
        @test isapprox(r.r_melt, 12375.0; atol=1e-9)
        @test isapprox(r.r_melt_ref, 13130.0; atol=1e-9)
        @test isapprox(r.r_melt - r.r_melt_ref, -755.0; atol=1e-9)
        @test isapprox(r.r_melt_error_cells, -1.01; atol=0.005)
        @test isapprox(r.tmax_rel_l2, 0.0237; atol=1e-4)
        @test r.energy_lse_residual <= 1.0e-9
        @test isapprox(bench.metrics["enthalpy_residual_1d"], 0.0061; atol=5e-5)
        # The pinned 1-D references of the suite record come from the same reference
        v129 = MSV.reference_values(bench.reference, 129)
        pinned(id) = MSV.suite_reference_value(id)
        @test isapprox(v129.centre_tmax, pinned("centre_tmax_nx129"); rtol=1e-12)
        @test isapprox(v129.melt_radius, pinned("melt_radius_nx129"); rtol=1e-12)
    end

    @testset "CI record: reproducibility, tampering and export round trip" begin
        ok, errs = MSV.check_validation_reproducibility(bench, MS_RECORD)
        @test ok
        @test isempty(errs)
        @test run_melting_shell_validation_check(MS_RECORD; bench=bench)
        @test run_validation_cli(["--check", MS_RECORD]; melting_shell_bench=bench) == 0
        @test last(VALIDATION_CASES).id == "melting_shell"
        raw = TOML.parsefile(MS_RECORD)
        @test [String(m["id"]) for m in raw["metrics"]] == collect(MSV.EXPECTED_METRIC_IDS)
        @test [Int(m["nx"]) for m in raw["metrics"]] == [65, 65, 65, 65, 0]
        mktempdir() do dir
            i_c = findfirst(m -> m["id"] == "centre_tmax_nx65", raw["metrics"])
            wide = deepcopy(raw)
            wide["metrics"][i_c]["atol"] = 40.0
            ok_w, errs_w = MSV.check_validation_reproducibility(
                bench, _ms_write_toml(joinpath(dir, "wide.toml"), wide)
            )
            @test !ok_w
            @test any(occursin("differs from the limit of the code", e) for e in errs_w)
            moved = deepcopy(raw)
            moved["metrics"][i_c]["value"] += 5.0
            moved["metadata"]["julia_version"] = string(VERSION)
            ok_m, errs_m = MSV.check_validation_reproducibility(
                bench, _ms_write_toml(joinpath(dir, "moved.toml"), moved)
            )
            @test !ok_m
            @test any(occursin("differs from recorded value", e) for e in errs_m)
            # A measured centre 25 K low fails the limit and the regression band
            low_metrics = copy(bench.metrics)
            low_metrics["centre_tmax_nx65"] -= 25.0
            low = (; bench..., metrics=low_metrics)
            ok_l, errs_l = MSV.check_validation_reproducibility(low, MS_RECORD)
            @test !ok_l
            @test any(occursin("fails its limit", e) for e in errs_l)
            @test any(occursin("physics regression band", e) for e in errs_l)
            # A centre 19.5 K high (a run at 0.75 L, measured once) passes the limit,
            # but the regression band of 0.382 K catches it
            @test isapprox(MSV.SEED_BAND_CENTRE_TMAX, 0.382; rtol=1e-14)
            hot_metrics = copy(bench.metrics)
            hot_metrics["centre_tmax_nx65"] += 19.5
            hot = (; bench..., metrics=hot_metrics)
            ok_h, errs_h = MSV.check_validation_reproducibility(hot, MS_RECORD)
            @test !ok_h
            @test !any(occursin("fails its limit", e) for e in errs_h)
            @test any(occursin("physics regression band", e) for e in errs_h)

            out = joinpath(dir, "exported.toml")
            commit = String(raw["metadata"]["src_commit"])
            export_cli(args, b) = MeltingShellExport.run_export_cli(
                args; src_commit=commit, bench=b, check_clean_tree=false
            )
            @test _ms_capture(() -> export_cli([out], bench)).value == 0
            @test first(MSV.check_validation_reproducibility(bench, out))
            @test occursin("regression band", read(out, String))
            @test _ms_capture(() -> export_cli(["--bad", out], bench)).value == 2
            refused = _ms_capture(() -> export_cli([out], low))
            @test refused.value == 1
            @test occursin("accept-physics-change", refused.err)
            @test_throws DomainError MeltingShellExport.run_export_cli(
                [out]; src_commit="abc", bench=bench, check_clean_tree=false
            )
        end
        @test !first(MSV.check_validation_reproducibility(bench, MS_RECORD * ".x"))
        @test_throws DomainError MSV.check_validation_reproducibility(bench, "")
    end
end
