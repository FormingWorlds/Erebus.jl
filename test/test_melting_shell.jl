using Test
using Erebus
isdefined(@__MODULE__, :MeltingShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_melting_shell.jl"))
const MSV = MeltingShellValidation

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
        for L in (4.0e5, 6.0e5)
            cfg = MSV.make_melting_shell_config(;
                Nx=17, phi_init=0.25, t_form_yr=1.3e6, n_steps=1, dt_yr=5.0e4, L_melt=L
            )
            @test cfg.melting.active
            @test cfg.melting.rho_melt == rho_s
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
        # Ste != 1 at L 6e5: latent term 4.95e6 differs from rho cp 3.3e6
        @test 3300.0 * 6.0e5 / 400.0 != 3.3e6
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
        @test hr[1] == hr[2]
    end

    @testset "Reference heat capacity matches the Erebus law" begin
        for L in (4.0e5, 6.0e5),
            phi in (0.0, 0.4),
            T in (1000.0, 1400.0, 1401.0, 1799.0, 1800.0, 2000.0)

            src = Erebus.rhocp_apparent_silicate(
                T, 0.0, 3.3e6, 3300.0, 1; T_sol=1400.0, T_liq=1800.0, L_melt=L
            )
            @test isapprox(MSV.rhocp_app_ref(T, phi; L_melt=L), (1 - phi) * src; rtol=1e-14)
        end
        @test MSV.melt_fraction_ref(1300.0) == 0.0
        @test MSV.melt_fraction_ref(1500.0) == 0.25
        @test MSV.melt_fraction_ref(1900.0) == 1.0
    end

    @testset "1-D reference: L = 0 equals the run without melting" begin
        kw = (;
            R=2.0e4, phi_init=0.25, t_form_yr=1.3e6, t_end_yr=2.0e6, dt_yr=5.0e4, Nr=200
        )
        a = MSV.PSV.solve_porous_sphere(; kw...)
        b = MSV.PSV.solve_porous_sphere(;
            kw..., melt=(; L_melt=0.0, T_sol=1400.0, T_liq=1800.0)
        )
        @test a.T == b.T && a.T_max == b.T_max && a.phi == b.phi && a.d_shell == b.d_shell
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
        @test r[1].T_peak > 1450.0
        e = [x.enthalpy_residual for x in r]
        @test e[1] > 1e-6
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
        @test r[1].F_max[1] == MSV.melt_fraction_ref(r[1].T_max[1])
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
        @test f(r -> r < 7500.0 ? 1600.0 : 1000.0) == 7500.0
        @test f(r -> 1000.0) == 0.0
        @test f(r -> 1600.0) == R_p
        # A hot ring outside a cold bin does not extend the radius
        @test f(r -> (r < 5000.0 || 6000.0 <= r < 7000.0) ? 1600.0 : 1000.0) == 5000.0
        # Exactly at the solidus a bin is not molten
        @test f(r -> r < 7500.0 ? 1600.0 : 1400.0) == 7500.0
        @test f(r -> r < 7500.0 ? 1400.0 : 1000.0) == 0.0
        # Empty bins are skipped
        keep = .!(3000.0 .<= rr .< 3500.0)
        @test MSV.melt_radius(
            xm[keep],
            ym[keep],
            [r < 7500.0 ? 1600.0 : 1000.0 for r in rr[keep]];
            xc=0.0,
            yc=0.0,
            R_p=R_p,
            dr=dr,
        ) == 7500.0
        # Median over sectors: 8 sectors to 7.5 km, 8 sectors to 8.5 km
        ang = mod.(atan.(ym, xm), 2pi)
        T2 = [
            (a < pi ? r < 7500.0 : r < 8500.0) ? 1600.0 : 1000.0 for (r, a) in zip(rr, ang)
        ]
        @test MSV.melt_radius(xm, ym, T2; xc=0.0, yc=0.0, R_p=R_p, dr=dr) == 8000.0
        # R_p / dr not an integer: the outer bin ends at R_p; markers at r >= R_p are ignored
        R2 = 2.02e4
        pts2 = lattice(R2 + 400.0)
        x2, y2 = first.(pts2), last.(pts2)
        r2 = hypot.(x2, y2)
        T3 = [r < R2 ? 1600.0 : 1000.0 for r in r2]
        @test MSV.melt_radius(x2, y2, T3; xc=0.0, yc=0.0, R_p=R2, dr=dr) == R2
        # Offset centre
        @test MSV.melt_radius(
            xm .+ 1.0e4,
            ym .- 3.0e3,
            [r < 7500.0 ? 1600.0 : 1000.0 for r in rr];
            xc=1.0e4,
            yc=-3.0e3,
            R_p=R_p,
            dr=dr,
        ) == 7500.0
        @test_throws ArgumentError MSV.melt_radius(
            [0.0], [0.0, 1.0], [1.0]; xc=0.0, yc=0.0, R_p=R_p, dr=dr
        )
        @test_throws DomainError MSV.melt_radius(
            [0.0], [0.0], [1.0]; xc=0.0, yc=0.0, R_p=R_p, dr=0.0
        )
    end
end
