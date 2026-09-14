using Test
using Erebus
using StaticArrays
using LinearAlgebra

include("test_helpers.jl")

@testset "Magma Ocean Multi-Component Degassing Physics" begin
    # Reference planet parameters (lunar-sized / Vesta-sized planetesimal)
    R_p = 50_000.0 # 50 km radius
    M_p = (4.0 / 3.0) * π * (R_p^3) * 3300.0 # ~1.73e18 kg
    g_surf = Erebus.G_GRAV * M_p / (R_p^2)
    T_mo = 1800.0 # Liquidus / super-liquidus temperature [K]
    M_melt = 0.8 * M_p # 80% magma ocean by mass

    # Test volatile inventories (representative carbonaceous chondrite volatile content)
    # H: 5000 ppm (as H2O equiv ~4.5 wt%)
    # C: 2000 ppm
    # N: 200 ppm
    # S: 10000 ppm (1 wt%)
    M_tot_H = M_p * 5.0e-3
    M_tot_C = M_p * 2.0e-3
    M_tot_N = M_p * 2.0e-4
    M_tot_S = M_p * 1.0e-2

    # -------------------------------------------------------------------------
    # 1. Equilibrium Magma Ocean Partitioning & Mass Conservation
    # -------------------------------------------------------------------------
    @testset "Equilibrium Volatile Partitioning & Mass Conservation" begin
        # Solve partitioning at neutral IW (delta_IW = 0.0)
        sol = solve_magma_ocean_volatile_partitioning(
            M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )

        @test sol isa NamedTuple
        @test sol.P_surf > 0.0
        @test isfinite(sol.P_surf)

        # Dalton's law: sum of partial pressures equals total surface pressure
        P_sum = sum(values(sol.p_i))
        @test isapprox(P_sum, sol.P_surf; rtol=1e-10)

        # Strict elemental mass conservation across melt and atmosphere:
        # 1. Hydrogen conservation (elemental H mass in melt + atmosphere)
        # H2O is 2.01588 / 18.01528 H by mass; H2 is 1.0 H by mass; CH4 is 4.032 / 16.042 H; NH3 is 3.024 / 17.031 H; H2S is 2.016 / 34.08 H
        M_atm_H = (
            get(sol.M_atm_i, :H2, 0.0) +
            get(sol.M_atm_i, :H2O, 0.0) * (2.01588 / 18.01528) +
            get(sol.M_atm_i, :CH4, 0.0) * (4.03176 / 16.04246) +
            get(sol.M_atm_i, :NH3, 0.0) * (3.02382 / 17.03052) +
            get(sol.M_atm_i, :H2S, 0.0) * (2.01588 / 34.08088)
        )
        M_melt_H = sol.M_melt_H
        @test isapprox(M_melt_H + M_atm_H, M_tot_H; rtol=1e-6)

        # 2. Carbon conservation
        M_atm_C = (
            get(sol.M_atm_i, :CO, 0.0) * (12.011 / 28.0101) +
            get(sol.M_atm_i, :CO2, 0.0) * (12.011 / 44.0095) +
            get(sol.M_atm_i, :CH4, 0.0) * (12.011 / 16.04246)
        )
        M_melt_C = sol.M_melt_C
        @test isapprox(M_melt_C + M_atm_C, M_tot_C; rtol=1e-6)

        # 3. Nitrogen conservation
        M_atm_N = (
            get(sol.M_atm_i, :N2, 0.0) + get(sol.M_atm_i, :NH3, 0.0) * (14.007 / 17.03052)
        )
        M_melt_N = sol.M_melt_N
        @test isapprox(M_melt_N + M_atm_N, M_tot_N; rtol=1e-6)

        # 4. Sulfur conservation
        M_atm_S = (
            get(sol.M_atm_i, :H2S, 0.0) * (32.06 / 34.08088) +
            get(sol.M_atm_i, :SO2, 0.0) * (32.06 / 64.066) +
            get(sol.M_atm_i, :S2, 0.0)
        )
        M_melt_S = sol.M_melt_S
        @test isapprox(M_melt_S + M_atm_S, M_tot_S; rtol=1e-6)

        # Monotonicity with total volatile inventory
        sol_2x = solve_magma_ocean_volatile_partitioning(
            M_melt,
            2.0 * M_tot_H,
            2.0 * M_tot_C,
            2.0 * M_tot_N,
            2.0 * M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )
        @test sol_2x.P_surf > sol.P_surf
        @test sol_2x.M_atm_tot > sol.M_atm_tot
        @test sol_2x.M_melt_tot > sol.M_melt_tot

        # Monotonicity with melt mass (more melt volume dissolves more volatiles, lowering surface pressure)
        sol_more_melt = solve_magma_ocean_volatile_partitioning(
            2.0 * M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )
        @test sol_more_melt.P_surf < sol.P_surf
        @test sol_more_melt.M_melt_tot > sol.M_melt_tot
    end

    # -------------------------------------------------------------------------
    # 2. Redox State Control on Degassed Volatile Speciation
    # -------------------------------------------------------------------------
    @testset "Redox Sensitivity on Outgassed Speciation" begin
        # Reduced case: delta_IW = -3.0 (enstatite/ordinary chondrite conditions)
        sol_red = solve_magma_ocean_volatile_partitioning(
            M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            -3.0;
            carbon_active=true,
            sulfur_active=true,
        )

        # Oxidized case: delta_IW = +2.0 (oxidized carbonaceous chondrite / terrestrial)
        sol_ox = solve_magma_ocean_volatile_partitioning(
            M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            2.0;
            carbon_active=true,
            sulfur_active=true,
        )

        # Reduced atmospheres must have higher CO/CO2 and H2/H2O ratios
        co_co2_red = sol_red.p_i[:CO] / max(sol_red.p_i[:CO2], 1e-30)
        co_co2_ox = sol_ox.p_i[:CO] / max(sol_ox.p_i[:CO2], 1e-30)
        @test co_co2_red > co_co2_ox

        h2_h2o_red = sol_red.p_i[:H2] / max(sol_red.p_i[:H2O], 1e-30)
        h2_h2o_ox = sol_ox.p_i[:H2] / max(sol_ox.p_i[:H2O], 1e-30)
        @test h2_h2o_red > h2_h2o_ox

        # Under reducing conditions, nitrogen dissolves chemically as nitride in melt
        # Hence, dissolved N in melt must be higher under reducing conditions
        @test sol_red.M_melt_N > sol_ox.M_melt_N
    end

    # -------------------------------------------------------------------------
    # 3. Asymptotic Limits & Boundary Cases
    # -------------------------------------------------------------------------
    @testset "Asymptotic Limits & Boundary Conditions" begin
        # Zero melt mass -> All volatiles in atmosphere
        sol_no_melt = solve_magma_ocean_volatile_partitioning(
            0.0,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )
        @test sol_no_melt.M_melt_tot ≈ 0.0 atol=1e-12
        @test isapprox(
            sol_no_melt.M_atm_tot, M_tot_H + M_tot_C + M_tot_N + M_tot_S; rtol=1e-4
        )

        # Zero volatile inventory -> Zero surface pressure
        sol_dry = solve_magma_ocean_volatile_partitioning(
            M_melt, 0.0, 0.0, 0.0, 0.0, R_p, g_surf, T_mo, 0.0
        )
        @test iszero(sol_dry.P_surf)
        @test iszero(sol_dry.M_atm_tot)

        # Error handling on invalid inputs
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            -1.0, M_tot_H, M_tot_C, M_tot_N, M_tot_S, R_p, g_surf, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, -1.0, M_tot_C, M_tot_N, M_tot_S, R_p, g_surf, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, M_tot_H, M_tot_C, M_tot_N, M_tot_S, 0.0, g_surf, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, M_tot_H, M_tot_C, M_tot_N, M_tot_S, R_p, 0.0, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, M_tot_H, M_tot_C, M_tot_N, M_tot_S, R_p, g_surf, 0.0, 0.0
        )
    end

    # -------------------------------------------------------------------------
    # 4. Marker-Scale Magma Ocean Degassing & Crystallization Exsolution
    # -------------------------------------------------------------------------
    @testset "Lagrangian Marker Magma Ocean Degassing" begin
        # Setup mock marker arrays
        N_markers = 100
        xm = fill(0.0, N_markers)
        ym = fill(0.95 * R_p, N_markers) # Near-surface markers
        tm = fill(2, N_markers) # Rock markers
        tkm = fill(1800.0, N_markers) # Molten
        Fm = fill(0.8, N_markers) # 80% melt fraction
        Fm_old = fill(0.8, N_markers)

        # High volatile concentration (supersaturated relative to low surface pressure)
        XH2Om = fill(0.02, N_markers) # 2 wt% H2O
        XCm = fill(1000.0e-6, N_markers) # 1000 ppm C
        XNm = fill(100.0e-6, N_markers) # 100 ppm N
        XSm = fill(2000.0e-6, N_markers) # 2000 ppm S
        rhosolidm_val = 3000.0
        marker_vol = ((2.0 * R_p)^2) / 10000.0 # Marker area in 2D

        cfg_degas = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            F_melt_threshold=0.40,
            degas_depth_fraction=0.90,
            crystallization_degassing=true,
            efficiency=1.0,
        )

        dt = 1000.0 # 1000 s
        P_surf = 1.0e5 # 1 bar

        # Decompression degassing pass
        rates = degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            Fm_old,
            XH2Om,
            XCm,
            XNm,
            XSm,
            N_markers,
            dt,
            P_surf,
            R_p,
            cfg_degas;
            rho_solid=rhosolidm_val,
            marker_volume=marker_vol,
            delta_IW=0.0,
        )

        @test rates isa Dict{Symbol,Float64}
        @test haskey(rates, :H2O)
        @test rates[:H2O] > 0.0

        # Marker volatiles must have decreased due to degassing
        @test all(XH2Om .< 0.02)

        # Crystallization exsolution pass: melt fraction decreases from 0.8 to 0.2
        Fm_crystallizing = fill(0.2, N_markers)
        Fm_prev = fill(0.8, N_markers)
        XH2O_pre_cryst = copy(XH2Om)

        rates_cryst = degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm_crystallizing,
            Fm_prev,
            XH2Om,
            XCm,
            XNm,
            XSm,
            N_markers,
            dt,
            P_surf,
            R_p,
            cfg_degas;
            rho_solid=rhosolidm_val,
            marker_volume=marker_vol,
            delta_IW=0.0,
        )

        @test rates_cryst[:H2O] > 0.0
        @test all(XH2Om .<= XH2O_pre_cryst)
    end

    # -------------------------------------------------------------------------
    # 5. Magma Ocean Degassing Integration in simulation_loop
    # -------------------------------------------------------------------------
    @testset "Magma Ocean Degassing Runtime Integration in simulation_loop" begin
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg_base = load_config(quick_toml)

        for mode in (:dynamic_flux, :equilibrium)
            mktempdir() do output_dir
                cfg_run = SimulationConfig(;
                    grid=cfg_base.grid,
                    geometry=cfg_base.geometry,
                    time=TimeConfig(
                        dt_initial=cfg_base.time.dt_initial,
                        dt_longest=cfg_base.time.dt_longest,
                        dtcoefdn=cfg_base.time.dtcoefdn,
                        dtcoefup=cfg_base.time.dtcoefup,
                        dtstep=cfg_base.time.dtstep,
                        dxymax=cfg_base.time.dxymax,
                        vpratio=cfg_base.time.vpratio,
                        DTmax=cfg_base.time.DTmax,
                        start_time=cfg_base.time.start_time,
                        endtime=cfg_base.time.endtime,
                        start_step=1,
                        n_steps=2,
                    ),
                    solver=cfg_base.solver,
                    poroelasticity=cfg_base.poroelasticity,
                    thermodynamics=cfg_base.thermodynamics,
                    reaction=cfg_base.reaction,
                    materials=cfg_base.materials,
                    output=OutputConfig(output_dir=output_dir, savematstep=2),
                    disk=cfg_base.disk,
                    melting=cfg_base.melting,
                    venting=VentingConfig(active=false),
                    magma_degassing=MagmaOceanDegassingConfig(active=true, mode=mode),
                    atmosphere=AtmosphereConfig(active=true, mode=:guillot),
                )

                Erebus.simulation_loop(cfg_run; output_path=output_dir)

                files = readdir(output_dir)
                @test "output_00000.jld2" in files
                @test "output_00002.jld2" in files

                data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
                @test data2["timestep"] == 2
                @test haskey(data2, "M_planet_val")
                @test data2["M_planet_val"] > 0.0
                @test haskey(data2, "atm_P_surf")
            end
        end
    end
end
