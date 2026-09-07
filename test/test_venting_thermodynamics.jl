using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using StaticArrays
using TOML

@testset "Venting Thermodynamics & Disk Dispersal" begin
    @testset "DiskConfig Dispersal Fields & Bounds Validation" begin
        cfg_def = DiskConfig()
        @test isapprox(cfg_def.t_dispersal_myr, 3.0; rtol=1e-12)
        @test isapprox(cfg_def.dt_dispersal_myr, 0.1; rtol=1e-12)
        @test isapprox(cfg_def.p_amb_disk, 10.0; rtol=1e-12)
        @test isapprox(cfg_def.p_amb_space, 1.0e-4; rtol=1e-12)
        @test isapprox(cfg_def.albedo, 0.06; rtol=1e-12)
        @test isnan(cfg_def.t_eq_custom)

        # Inclusion in SimulationConfig
        sim_cfg = default_config()
        @test sim_cfg.disk isa DiskConfig
        @test validate_config(sim_cfg) === nothing

        # Validation: non-positive dispersal time
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; t_dispersal_myr=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; t_dispersal_myr=-1.0))
        )

        # Validation: non-positive transition duration
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; dt_dispersal_myr=0.0))
        )

        # Validation: non-positive disk ambient pressure
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; p_amb_disk=0.0))
        )

        # Validation: negative space vacuum floor
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; p_amb_space=-1.0e-6))
        )

        # Validation: unphysical albedo
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; albedo=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; albedo=1.0))
        )

        # Validation: unphysical custom equilibrium temperature
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; t_eq_custom=-50.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; disk=DiskConfig(; t_eq_custom=0.0))
        )
    end

    @testset "VentingConfig Schema & TOML Serialization" begin
        v_def = VentingConfig()
        @test v_def.active == false
        @test v_def.mode === :darcy_sink
        @test isapprox(v_def.k_vent, 1.0e-11; rtol=1e-12)
        @test isapprox(v_def.conductance_factor, 1.0; rtol=1e-12)
        @test isapprox(v_def.L_sublimation, 2.83e6; rtol=1e-12)
        @test v_def.latent_cooling == true

        # Top-level inclusion
        sim_cfg = default_config()
        @test sim_cfg.venting isa VentingConfig
        @test validate_config(sim_cfg) === nothing

        # Validation: invalid mode
        @test_throws ArgumentError validate_config(
            SimulationConfig(; venting=VentingConfig(; mode=:explosive_breakup))
        )

        # Validation: non-positive permeability
        @test_throws ArgumentError validate_config(
            SimulationConfig(; venting=VentingConfig(; k_vent=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; venting=VentingConfig(; k_vent=-1.0e-12))
        )

        # Validation: non-positive conductance factor
        @test_throws ArgumentError validate_config(
            SimulationConfig(; venting=VentingConfig(; conductance_factor=0.0))
        )

        # Validation: non-positive latent heat of sublimation
        @test_throws ArgumentError validate_config(
            SimulationConfig(; venting=VentingConfig(; L_sublimation=0.0))
        )

        # TOML round-trip serialization
        cfg_custom = SimulationConfig(;
            venting=VentingConfig(;
                active=true,
                mode=:darcy_sink,
                k_vent=5.0e-11,
                conductance_factor=2.0,
                L_sublimation=2.8e6,
                latent_cooling=false,
            ),
        )
        toml_str = save_config(cfg_custom)
        parsed_cfg = load_config(toml_str)
        @test parsed_cfg.venting.active == true
        @test parsed_cfg.venting.mode === :darcy_sink
        @test isapprox(parsed_cfg.venting.k_vent, 5.0e-11; rtol=1e-12)
        @test isapprox(parsed_cfg.venting.conductance_factor, 2.0; rtol=1e-12)
        @test isapprox(parsed_cfg.venting.L_sublimation, 2.8e6; rtol=1e-12)
        @test parsed_cfg.venting.latent_cooling == false
    end

    @testset "Disk Dispersal Transition Weight" begin
        # Halfway point at dispersal epoch
        w_half = compute_disk_dispersal_weight(
            3.0e6 * 365.25 * 86400.0; t_dispersal_myr=3.0, dt_dispersal_myr=0.1
        )
        @test isapprox(w_half, 0.5; atol=1e-6)

        # Well before dispersal
        w_early = compute_disk_dispersal_weight(
            2.0e6 * 365.25 * 86400.0; t_dispersal_myr=3.0, dt_dispersal_myr=0.1
        )
        @test isapprox(w_early, 0.0; atol=1e-4)

        # Well after dispersal
        w_late = compute_disk_dispersal_weight(
            4.0e6 * 365.25 * 86400.0; t_dispersal_myr=3.0, dt_dispersal_myr=0.1
        )
        @test isapprox(w_late, 1.0; atol=1e-4)

        # Monotonicity
        t_seq = range(2.5, 3.5; length=11) .* (1.0e6 * 365.25 * 86400.0)
        w_seq = [
            compute_disk_dispersal_weight(t; t_dispersal_myr=3.0, dt_dispersal_myr=0.1) for
            t in t_seq
        ]
        @test issorted(w_seq)
        @test all(w -> 0.0 <= w <= 1.0, w_seq)

        # Error guards
        @test_throws DomainError compute_disk_dispersal_weight(
            1.0; t_dispersal_myr=3.0, dt_dispersal_myr=0.0
        )
        @test_throws DomainError compute_disk_dispersal_weight(
            1.0; t_dispersal_myr=3.0, dt_dispersal_myr=-0.1
        )
    end

    @testset "Solar Radiation Equilibrium Temperature" begin
        # 2.7 AU (asteroid belt): Teq ~ 166.8 K
        T_eq_27 = compute_solar_equilibrium_temperature(2.7; albedo=0.06)
        @test 160.0 <= T_eq_27 <= 172.0
        @test isapprox(T_eq_27, 166.8; atol=2.0)

        # 1.0 AU (Earth orbit): Teq ~ 274 K (for A=0.06)
        T_eq_10 = compute_solar_equilibrium_temperature(1.0; albedo=0.06)
        @test 270.0 <= T_eq_10 <= 280.0

        # Radial scaling: Teq proportional to d^(-0.5)
        ratio = T_eq_10 / T_eq_27
        expected_ratio = sqrt(2.7 / 1.0)
        @test isapprox(ratio, expected_ratio; rtol=1e-4)

        # Error guards
        @test_throws DomainError compute_solar_equilibrium_temperature(0.0)
        @test_throws DomainError compute_solar_equilibrium_temperature(-1.0)
        @test_throws DomainError compute_solar_equilibrium_temperature(1.0; albedo=-0.05)
        @test_throws DomainError compute_solar_equilibrium_temperature(1.0; albedo=1.0)
    end

    @testset "Ambient Conditions Evolution Across Dispersal" begin
        # Default: dispersal_active = false preserves constant nebular ambient
        disk_cfg_nodisp = DiskConfig(;
            enabled=true, model=:fixed, t_ambient=170.0, p_amb_disk=10.0
        )
        yr_sec = 365.25 * 86400.0
        T_nodisp, P_nodisp, w_nodisp = compute_ambient_conditions(
            5.0e6 * yr_sec, disk_cfg_nodisp
        )
        @test isapprox(T_nodisp, 170.0; atol=1e-12)
        @test isapprox(P_nodisp, 10.0; atol=1e-12)
        @test iszero(w_nodisp)

        disk_cfg = DiskConfig(;
            enabled=true,
            model=:fixed,
            t_ambient=170.0,
            orbital_distance_au=2.7,
            t_dispersal_myr=3.0,
            dt_dispersal_myr=0.1,
            p_amb_disk=10.0,
            p_amb_space=1.0e-4,
            albedo=0.06,
            dispersal_active=true,
        )

        # Before dispersal: nebular values
        T_pre, P_pre, w_pre = compute_ambient_conditions(1.0e6 * yr_sec, disk_cfg)
        @test isapprox(T_pre, 170.0; atol=1e-2)
        @test isapprox(P_pre, 10.0; atol=1e-2)
        @test isapprox(w_pre, 0.0; atol=1e-3)

        # At dispersal epoch: average values
        T_mid, P_mid, w_mid = compute_ambient_conditions(3.0e6 * yr_sec, disk_cfg)
        @test isapprox(w_mid, 0.5; atol=1e-4)
        @test isapprox(P_mid, 5.00005; atol=1e-2)
        T_eq_expected = compute_solar_equilibrium_temperature(2.7; albedo=0.06)
        @test isapprox(T_mid, 0.5 * (170.0 + T_eq_expected); atol=1e-2)

        # Post dispersal: vacuum space values
        T_post, P_post, w_post = compute_ambient_conditions(5.0e6 * yr_sec, disk_cfg)
        @test isapprox(w_post, 1.0; atol=1e-3)
        @test isapprox(P_post, 1.0e-4; atol=1e-5)
        @test isapprox(T_post, T_eq_expected; atol=1e-2)

        # Custom equilibrium temperature override
        disk_cfg_custom = DiskConfig(;
            enabled=true,
            model=:fixed,
            t_ambient=200.0,
            t_dispersal_myr=2.0,
            dt_dispersal_myr=0.05,
            t_eq_custom=150.0,
            dispersal_active=true,
        )
        T_c, _, _ = compute_ambient_conditions(4.0e6 * yr_sec, disk_cfg_custom)
        @test isapprox(T_c, 150.0; atol=1e-2)
    end

    @testset "Clausius-Clapeyron Ice Sublimation Vapor Pressure" begin
        # Triple point: T = 273.16 K -> P = 611.66 Pa
        P_triple = compute_ice_vapor_pressure(273.16)
        @test isapprox(P_triple, 611.66; rtol=1e-3)
        @test isapprox(compute_water_vapor_pressure(273.16), 611.66; rtol=1e-3)

        # Cold lid temperatures
        P_200 = compute_ice_vapor_pressure(200.0)
        @test 0.10 <= P_200 <= 0.25
        @test isapprox(P_200, 0.165; rtol=0.05)

        P_150 = compute_ice_vapor_pressure(150.0)
        @test 1.0e-6 <= P_150 <= 2.0e-5
        @test isapprox(P_150, 5.9e-6; rtol=0.1)

        # Liquid regime and boiling point: T = 373.15 K -> P ~ 101.3 kPa (1 atm)
        P_boil = compute_water_vapor_pressure(373.15)
        @test isapprox(P_boil, 101325.0; rtol=0.01)

        # High-temperature hydrothermal regime: T = 500 K -> P ~ 2.46 MPa
        P_500 = compute_water_vapor_pressure(500.0)
        @test isapprox(P_500, 2.46e6; rtol=0.05)

        # Supercritical regime clamp
        P_crit = compute_water_vapor_pressure(700.0)
        @test isapprox(P_crit, 22.064e6; rtol=1e-12)

        # Monotonicity across all physical regimes
        temps = [140.0, 180.0, 220.0, 273.16, 300.0, 373.15, 450.0, 550.0, 647.096]
        pressures = [compute_water_vapor_pressure(T) for T in temps]
        @test issorted(pressures)
        @test all(p -> p > 0.0, pressures)

        # Domain bounds and parameter guards
        @test_throws DomainError compute_water_vapor_pressure(0.0)
        @test_throws DomainError compute_water_vapor_pressure(-10.0)
        @test_throws DomainError compute_water_vapor_pressure(200.0; P0=0.0)
        @test_throws DomainError compute_water_vapor_pressure(200.0; T0=-1.0)
        @test_throws DomainError compute_water_vapor_pressure(200.0; Rv=0.0)
    end

    @testset "Venting Pressure with Cold-Trap Limit" begin
        # Cold-trap regime: ambient nebular gas exceeds sublimation vapor pressure
        # T_surf = 170 K -> P_sat ~ 1.6e-3 Pa; P_amb = 10 Pa
        P_vent_trap = compute_venting_pressure(170.0, 10.0)
        @test isapprox(P_vent_trap, 10.0; rtol=1e-12)
        @test P_vent_trap >= compute_ice_vapor_pressure(170.0)

        # Sublimation-dominated regime: warm rock in thin gas / vacuum
        # T_surf = 250 K -> P_sat ~ 76 Pa; P_amb = 1.0 Pa
        P_vent_sub = compute_venting_pressure(250.0, 1.0)
        P_sat_250 = compute_ice_vapor_pressure(250.0)
        @test isapprox(P_vent_sub, P_sat_250; rtol=1e-12)
        @test P_vent_sub > 1.0

        # Vacuum space limit
        P_vent_vac = compute_venting_pressure(200.0, 1.0e-4)
        P_sat_200 = compute_ice_vapor_pressure(200.0)
        @test isapprox(P_vent_vac, P_sat_200; rtol=1e-12)
    end
end
