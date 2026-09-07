using Test
using Erebus
using Erebus.Config
using Erebus.Physics

@testset "Multi-Species Volatile Solubility & Nitrogen Chemistry" begin
    @testset "compute_iron_wustite_fO2 Invariants" begin
        # Standard temperatures
        T_1000K = 1000.0
        log_fO2_1000 = compute_iron_wustite_fO2(T_1000K)
        # Expected: 6.541 - 28164 / 1000 = 6.541 - 28.164 = -21.623 [log10 bar]
        @test isapprox(log_fO2_1000, -21.623; rtol=1e-5)

        # Monotonicity with temperature
        T_1200K = 1200.0
        log_fO2_1200 = compute_iron_wustite_fO2(T_1200K)
        @test log_fO2_1200 > log_fO2_1000

        # With delta_IW offset
        log_fO2_reduced = compute_iron_wustite_fO2(T_1000K; delta_IW=-2.0)
        @test isapprox(log_fO2_reduced, log_fO2_1000 - 2.0; rtol=1e-5)

        # DomainError guards on non-physical temperature
        @test_throws DomainError compute_iron_wustite_fO2(0.0)
        @test_throws DomainError compute_iron_wustite_fO2(-100.0)
        @test_throws DomainError compute_iron_wustite_fO2(NaN)
        @test_throws DomainError compute_iron_wustite_fO2(Inf)
    end

    @testset "compute_water_solubility_melt Invariants & Asymptotics" begin
        As = 0.40 # wt% / MPa^0.5

        # Zero or negative pressure: zero dissolved water
        @test iszero(compute_water_solubility_melt(0.0; As=As))
        @test iszero(compute_water_solubility_melt(-1.0e5; As=As))

        # 1 MPa (1.0e6 Pa): exactly As wt%
        w_1MPa = compute_water_solubility_melt(1.0e6; As=As)
        @test isapprox(w_1MPa, 0.40; rtol=1e-12)

        # 4 MPa (4.0e6 Pa): sqrt(4) * As = 2 * 0.4 = 0.80 wt%
        w_4MPa = compute_water_solubility_melt(4.0e6; As=As)
        @test isapprox(w_4MPa, 0.80; rtol=1e-12)

        # 100 MPa: sqrt(100) * As = 10 * 0.4 = 4.0 wt%
        w_100MPa = compute_water_solubility_melt(1.0e8; As=As)
        @test isapprox(w_100MPa, 4.0; rtol=1e-12)

        # Scaling law: quadrupling pressure doubles dissolved water
        @test isapprox(w_4MPa / w_1MPa, 2.0; rtol=1e-12)

        # DomainError guard on non-finite pressure
        @test_throws DomainError compute_water_solubility_melt(NaN)
        @test_throws DomainError compute_water_solubility_melt(Inf)
        @test_throws DomainError compute_water_solubility_melt(1.0e6; As=-0.1)
    end

    @testset "compute_nitrogen_solubility_melt Redox Scaling & Partitioning" begin
        # 10 MPa pore pressure (100 bar)
        P_10MPa = 10.0e6

        # Case 1: Neutral IW buffer (delta_IW = 0)
        sol_IW = compute_nitrogen_solubility_melt(P_10MPa, 0.0)
        @test sol_IW.total_ppm > 0.0
        @test sol_IW.physical_ppm > 0.0
        @test sol_IW.chemical_ppm > 0.0
        @test isapprox(
            sol_IW.total_ppm, sol_IW.physical_ppm + sol_IW.chemical_ppm; rtol=1e-12
        )

        # Keyword default call matches neutral IW
        sol_default = compute_nitrogen_solubility_melt(P_10MPa)
        @test isapprox(sol_default.total_ppm, sol_IW.total_ppm; rtol=1e-12)

        # Case 2: Reducing conditions (delta_IW = -2)
        # Chemical nitride solubility scales as fO2^(-3/4).
        # When delta_IW drops by 2 log units, fO2 drops by 10^2 = 100, so fO2^(-3/4) increases by 100^(0.75) = 10^1.5 ≈ 31.62
        sol_reduced = compute_nitrogen_solubility_melt(P_10MPa, -2.0)
        @test isapprox(sol_reduced.physical_ppm, sol_IW.physical_ppm; rtol=1e-10)
        ratio_nitride = sol_reduced.chemical_ppm / sol_IW.chemical_ppm
        expected_ratio = 10.0^(2.0 * 0.75) # 10^1.5 = 31.6227766
        @test isapprox(ratio_nitride, expected_ratio; rtol=1e-4)
        @test sol_reduced.total_ppm > sol_IW.total_ppm

        # Case 3: Oxidizing conditions (delta_IW = +4) -> physical dissolution dominates
        sol_oxidized = compute_nitrogen_solubility_melt(P_10MPa, 4.0)
        @test sol_oxidized.physical_ppm > sol_oxidized.chemical_ppm

        # Zero pressure limit
        sol_zero = compute_nitrogen_solubility_melt(0.0, 0.0)
        @test iszero(sol_zero.total_ppm)
        @test iszero(sol_zero.physical_ppm)
        @test iszero(sol_zero.chemical_ppm)

        # DomainError guards
        @test_throws DomainError compute_nitrogen_solubility_melt(NaN, 0.0)
        @test_throws DomainError compute_nitrogen_solubility_melt(P_10MPa, NaN)
        @test_throws DomainError compute_nitrogen_solubility_melt(P_10MPa, -400.0)
        @test_throws DomainError compute_nitrogen_solubility_melt(P_10MPa, 400.0)
        @test_throws DomainError compute_nitrogen_solubility_melt(P_10MPa; Kh=-0.1)
        @test_throws DomainError compute_nitrogen_solubility_melt(P_10MPa; C_nitride=-0.1)
    end

    @testset "compute_organic_nitrogen_yield Invariants" begin
        T_devol = 550.0
        dT = 50.0

        # Midpoint symmetry: at T = T_devol, yield is exactly 0.5
        @test isapprox(
            compute_organic_nitrogen_yield(T_devol; T_devol=T_devol, delta_T=dT),
            0.5;
            rtol=1e-12,
        )

        # Low temperature asymptote: T << T_devol -> 0
        y_cold = compute_organic_nitrogen_yield(200.0; T_devol=T_devol, delta_T=dT)
        @test y_cold < 2.0e-3
        @test y_cold >= 0.0

        # Deep cold asymptote (T = 100 K)
        y_deep_cold = compute_organic_nitrogen_yield(100.0; T_devol=T_devol, delta_T=dT)
        @test y_deep_cold < 1.0e-3
        @test y_deep_cold >= 0.0

        # High temperature asymptote: T >> T_devol -> 1
        y_hot = compute_organic_nitrogen_yield(900.0; T_devol=T_devol, delta_T=dT)
        @test isapprox(y_hot, 1.0; atol=1e-3)
        @test y_hot <= 1.0

        # Monotonicity
        @test compute_organic_nitrogen_yield(600.0) > compute_organic_nitrogen_yield(500.0)

        # DomainError guards
        @test_throws DomainError compute_organic_nitrogen_yield(0.0)
        @test_throws DomainError compute_organic_nitrogen_yield(-50.0)
        @test_throws DomainError compute_organic_nitrogen_yield(NaN)
        @test_throws DomainError compute_organic_nitrogen_yield(500.0; delta_T=0.0)
    end

    @testset "VolatilesConfig Schema & Bounds Validation" begin
        v_default = VolatilesConfig()
        @test !v_default.active
        @test isapprox(v_default.fO2_delta_IW, -1.0; rtol=1e-12)
        @test isapprox(v_default.water_solubility_coeff, 0.40; rtol=1e-12)
        @test isapprox(v_default.nitrogen_henry_coeff, 0.40; rtol=1e-12)

        # Valid config
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg_base = load_config(quick_toml)
        cfg_vol = SimulationConfig(
            grid=cfg_base.grid,
            geometry=cfg_base.geometry,
            time=cfg_base.time,
            solver=cfg_base.solver,
            poroelasticity=cfg_base.poroelasticity,
            thermodynamics=cfg_base.thermodynamics,
            reaction=cfg_base.reaction,
            materials=cfg_base.materials,
            output=cfg_base.output,
            disk=cfg_base.disk,
            melting=cfg_base.melting,
            venting=cfg_base.venting,
            volatiles=VolatilesConfig(
                active=true,
                fO2_delta_IW=-2.0,
                water_solubility_coeff=0.45,
                nitrogen_henry_coeff=0.50,
                nitrogen_nitride_capacity=2.0e-3,
                t_organic_devol=520.0,
                dt_organic_devol=40.0,
                organic_n_initial_ppm=600.0,
            ),
        )
        @test validate_config(cfg_vol) === nothing

        # Invalid bounds checks
        # Non-positive water_solubility_coeff
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(water_solubility_coeff=-0.1),
            ),
        )
        # Non-positive nitrogen_henry_coeff
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(nitrogen_henry_coeff=0.0),
            ),
        )
        # Non-positive nitrogen_nitride_capacity
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(nitrogen_nitride_capacity=-1.0e-4),
            ),
        )
        # Non-positive t_organic_devol
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(t_organic_devol=0.0),
            ),
        )
        # Non-positive dt_organic_devol
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(dt_organic_devol=-10.0),
            ),
        )
        # Negative organic_n_initial_ppm
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(organic_n_initial_ppm=-10.0),
            ),
        )
        # Out-of-bounds fO2_delta_IW
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(fO2_delta_IW=-100.0),
            ),
        )
    end

    @testset "VolatilesConfig Serialization Round-Trip" begin
        cfg_orig = SimulationConfig(
            volatiles=VolatilesConfig(
                active=true,
                fO2_delta_IW=-2.5,
                water_solubility_coeff=0.42,
                nitrogen_henry_coeff=0.48,
                nitrogen_nitride_capacity=2.5e-3,
                t_organic_devol=530.0,
                dt_organic_devol=45.0,
                organic_n_initial_ppm=650.0,
            ),
        )
        toml_str = save_config(cfg_orig)
        cfg_loaded = load_config(toml_str)
        @test cfg_loaded.volatiles.active == cfg_orig.volatiles.active
        @test isapprox(
            cfg_loaded.volatiles.fO2_delta_IW, cfg_orig.volatiles.fO2_delta_IW; rtol=1e-12
        )
        @test isapprox(
            cfg_loaded.volatiles.water_solubility_coeff,
            cfg_orig.volatiles.water_solubility_coeff;
            rtol=1e-12,
        )
        @test isapprox(
            cfg_loaded.volatiles.nitrogen_henry_coeff,
            cfg_orig.volatiles.nitrogen_henry_coeff;
            rtol=1e-12,
        )
        @test isapprox(
            cfg_loaded.volatiles.nitrogen_nitride_capacity,
            cfg_orig.volatiles.nitrogen_nitride_capacity;
            rtol=1e-12,
        )
        @test isapprox(
            cfg_loaded.volatiles.t_organic_devol,
            cfg_orig.volatiles.t_organic_devol;
            rtol=1e-12,
        )
        @test isapprox(
            cfg_loaded.volatiles.dt_organic_devol,
            cfg_orig.volatiles.dt_organic_devol;
            rtol=1e-12,
        )
        @test isapprox(
            cfg_loaded.volatiles.organic_n_initial_ppm,
            cfg_orig.volatiles.organic_n_initial_ppm;
            rtol=1e-12,
        )
    end
end
