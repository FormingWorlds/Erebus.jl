using Test
using Erebus
using Erebus.Physics
using Erebus.Particles
using TOML

@testset "Hydrothermal Subgrid Convection" begin
    @testset "HydrothermalConfig Defaults and Validation" begin
        cfg_def = HydrothermalConfig()
        @test cfg_def.active == false
        @test cfg_def.phi_start ≈ 0.30
        @test cfg_def.phi_end ≈ 0.70
        @test cfg_def.Ra_m_crit ≈ 4.0 * pi^2
        @test cfg_def.Ra_crit ≈ 1100.0
        @test cfg_def.c_porous ≈ 1.0
        @test cfg_def.c_free ≈ 0.088
        @test cfg_def.H_layer ≈ 10000.0
        @test cfg_def.dT_min ≈ 5.0
        @test cfg_def.k_floor ≈ 1.0e-3
        @test cfg_def.k_cutoff ≈ 1.0e6
        @test cfg_def.picard_damping ≈ 0.5
        @test cfg_def.resolution_weighting == true
        @test cfg_def.Pe_crit ≈ 2.0
        @test cfg_def.T_surface_ref ≈ 273.15
        @test cfg_def.gravity ≈ 0.5
        @test cfg_def.dynamic_gravity == true
        @test cfg_def.dynamic_layer_depth == true
        @test cfg_def.dynamic_peclet == true
        @test cfg_def.H_layer_min ≈ 100.0

        # TOML serialization and deserialization roundtrip
        sim_cfg = SimulationConfig(;
            hydrothermal=HydrothermalConfig(;
                active=true,
                phi_start=0.25,
                phi_end=0.65,
                picard_damping=0.75,
                c_free=0.10,
                dynamic_gravity=false,
                dynamic_layer_depth=false,
                dynamic_peclet=false,
                H_layer_min=250.0,
            ),
        )
        toml_str = save_config(sim_cfg)
        loaded_cfg = load_config(toml_str)
        @test loaded_cfg.hydrothermal.active == true
        @test loaded_cfg.hydrothermal.phi_start ≈ 0.25
        @test loaded_cfg.hydrothermal.phi_end ≈ 0.65
        @test loaded_cfg.hydrothermal.picard_damping ≈ 0.75
        @test loaded_cfg.hydrothermal.c_free ≈ 0.10
        @test loaded_cfg.hydrothermal.dynamic_gravity == false
        @test loaded_cfg.hydrothermal.dynamic_layer_depth == false
        @test loaded_cfg.hydrothermal.dynamic_peclet == false
        @test loaded_cfg.hydrothermal.H_layer_min ≈ 250.0

        # Unphysical parameter validation
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, phi_start=-0.1)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, phi_start=0.7, phi_end=0.5)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; hydrothermal=HydrothermalConfig(; active=true, phi_end=1.5))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, Ra_m_crit=-10.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; hydrothermal=HydrothermalConfig(; active=true, Ra_crit=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, c_porous=-1.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; hydrothermal=HydrothermalConfig(; active=true, c_free=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, H_layer=-500.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; hydrothermal=HydrothermalConfig(; active=true, dT_min=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; hydrothermal=HydrothermalConfig(; active=true, k_floor=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, k_floor=10.0, k_cutoff=5.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, picard_damping=0.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, picard_damping=1.2)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; hydrothermal=HydrothermalConfig(; active=true, Pe_crit=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(; active=true, H_layer_min=-10.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                hydrothermal=HydrothermalConfig(;
                    active=true, H_layer_min=20000.0, H_layer=10000.0
                ),
            ),
        )
    end

    @testset "Porous Rayleigh-Darcy Scaling" begin
        # Baseline fluid and physical properties
        rho_f = 1000.0
        cp_f = 4184.0
        g = 0.5
        alpha_f = 2.0e-4
        K = 1.0e-12
        dT = 50.0
        H = 10000.0
        mu_f = 1.0e-3
        k_cond = 2.5

        # Analytical value check
        # Ra_m = (1000^2 * 4184 * 0.5 * 2e-4 * 1e-12 * 50 * 10000) / (1e-3 * 2.5) = 83.68
        Ra_m = compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, dT, H, mu_f, k_cond
        )
        @test Ra_m ≈ 83.68 rtol = 1.0e-12

        # Inactive thermal contrast or permeability yields zero
        @test compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, 0.0, H, mu_f, k_cond
        ) ≈ 0.0 atol = 1.0e-14
        @test compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, -10.0, H, mu_f, k_cond
        ) ≈ 0.0 atol = 1.0e-14
        @test compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, 0.0, dT, H, mu_f, k_cond
        ) ≈ 0.0 atol = 1.0e-14

        # Linear scaling assertions
        Ra_m_2x_K = compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, 2.0 * K, dT, H, mu_f, k_cond
        )
        @test Ra_m_2x_K ≈ 2.0 * Ra_m rtol = 1.0e-12

        Ra_m_2x_dT = compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, 2.0 * dT, H, mu_f, k_cond
        )
        @test Ra_m_2x_dT ≈ 2.0 * Ra_m rtol = 1.0e-12

        Ra_m_half_mu = compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, dT, H, 0.5 * mu_f, k_cond
        )
        @test Ra_m_half_mu ≈ 2.0 * Ra_m rtol = 1.0e-12

        # Error guards
        @test_throws DomainError compute_porous_rayleigh_darcy(
            -rho_f, cp_f, g, alpha_f, K, dT, H, mu_f, k_cond
        )
        @test_throws DomainError compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, -K, dT, H, mu_f, k_cond
        )
        @test_throws DomainError compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, dT, H, 0.0, k_cond
        )
        @test_throws DomainError compute_porous_rayleigh_darcy(
            rho_f, cp_f, g, alpha_f, K, dT, H, mu_f, -k_cond
        )
        @test_throws DomainError compute_porous_rayleigh_darcy(
            NaN, cp_f, g, alpha_f, K, dT, H, mu_f, k_cond
        )
    end

    @testset "Free-Fluid Rayleigh Scaling" begin
        rho_f = 1000.0
        cp_f = 4184.0
        g = 0.5
        alpha_f = 2.0e-4
        dT = 20.0
        H = 5000.0
        mu_f = 1.0e-3
        k_f = 0.6

        # Analytical value check
        # Ra = (1000^2 * 4184 * 0.5 * 2e-4 * 20 * 5000^3) / (1e-3 * 0.6) = 1.7433333333333333e21
        Ra = compute_free_fluid_rayleigh(rho_f, cp_f, g, alpha_f, dT, H, mu_f, k_f)
        @test Ra ≈ 1.7433333333333333e21 rtol = 1.0e-12

        # Cubic layer thickness scaling
        Ra_2x_H = compute_free_fluid_rayleigh(
            rho_f, cp_f, g, alpha_f, dT, 2.0 * H, mu_f, k_f
        )
        @test Ra_2x_H ≈ 8.0 * Ra rtol = 1.0e-12

        # Sub-threshold temperature
        @test compute_free_fluid_rayleigh(rho_f, cp_f, g, alpha_f, 0.0, H, mu_f, k_f) ≈ 0.0 atol =
            1.0e-14
        @test compute_free_fluid_rayleigh(rho_f, cp_f, g, alpha_f, -5.0, H, mu_f, k_f) ≈ 0.0 atol =
            1.0e-14

        # Error guards
        @test_throws DomainError compute_free_fluid_rayleigh(
            rho_f, cp_f, g, alpha_f, dT, -H, mu_f, k_f
        )
        @test_throws DomainError compute_free_fluid_rayleigh(
            rho_f, cp_f, g, alpha_f, dT, H, 0.0, k_f
        )
        @test_throws DomainError compute_free_fluid_rayleigh(
            rho_f, cp_f, g, alpha_f, dT, H, mu_f, 0.0
        )
        @test_throws DomainError compute_free_fluid_rayleigh(
            rho_f, cp_f, g, alpha_f, dT, H, mu_f, Inf
        )
    end

    @testset "Hydrothermal Nusselt Number Blending" begin
        Ra_m_crit = 4.0 * pi^2
        Ra_crit = 1100.0

        # 1. Sub-critical conditions in both porous and free regimes
        @test compute_hydrothermal_nusselt(0.5 * Ra_m_crit, 0.5 * Ra_crit, 0.1) ≈ 1.0 atol =
            1.0e-14
        @test compute_hydrothermal_nusselt(0.5 * Ra_m_crit, 0.5 * Ra_crit, 0.9) ≈ 1.0 atol =
            1.0e-14
        @test compute_hydrothermal_nusselt(0.5 * Ra_m_crit, 0.5 * Ra_crit, 0.5) ≈ 1.0 atol =
            1.0e-14

        # 2. Pure porous regime (phi <= phi_start = 0.30)
        # Nu_porous = 1.0 + 1.0 * (10.0 * Ra_m_crit / Ra_m_crit - 1.0) = 10.0
        Nu_porous_10x = compute_hydrothermal_nusselt(10.0 * Ra_m_crit, 1.0e8, 0.20)
        @test Nu_porous_10x ≈ 10.0 rtol = 1.0e-12

        # 3. Pure free-fluid regime (phi >= phi_end = 0.70)
        # Ra = 1.0e6 => cbrt(1e6) = 100.0 => Nu_free = max(1.0, 0.088 * 100) = 8.8
        Nu_free_1e6 = compute_hydrothermal_nusselt(10.0 * Ra_m_crit, 1.0e6, 0.85)
        @test Nu_free_1e6 ≈ 8.8 rtol = 1.0e-12

        # 4. Smoothstep monotonicity and C^0 continuity across transition
        phi_vals = range(0.30, 0.70; length=50)
        Nu_vals = [
            compute_hydrothermal_nusselt(
                10.0 * Ra_m_crit, 1.0e6, p; phi_start=0.30, phi_end=0.70
            ) for p in phi_vals
        ]
        @test Nu_vals[1] ≈ 10.0 rtol = 1.0e-12
        @test Nu_vals[end] ≈ 8.8 rtol = 1.0e-12
        for k in 1:(length(Nu_vals) - 1)
            @test Nu_vals[k] >= Nu_vals[k + 1] - 1.0e-14
        end

        # Reversed case: Nu_porous < Nu_free
        Nu_rev_vals = [
            compute_hydrothermal_nusselt(
                2.0 * Ra_m_crit, 1.0e9, p; phi_start=0.30, phi_end=0.70
            ) for p in phi_vals
        ]
        @test Nu_rev_vals[1] ≈ 2.0 rtol = 1.0e-12
        @test Nu_rev_vals[end] ≈ 88.0 rtol = 1.0e-12
        for k in 1:(length(Nu_rev_vals) - 1)
            @test Nu_rev_vals[k] <= Nu_rev_vals[k + 1] + 1.0e-14
        end

        # 5. Nusselt lower bound is strictly 1.0
        @test compute_hydrothermal_nusselt(0.0, 0.0, 0.5) ≈ 1.0 atol = 1.0e-14

        # Error guards
        @test_throws DomainError compute_hydrothermal_nusselt(-1.0, 100.0, 0.5)
        @test_throws DomainError compute_hydrothermal_nusselt(10.0, -100.0, 0.5)
        @test_throws DomainError compute_hydrothermal_nusselt(10.0, 100.0, -0.1)
        @test_throws DomainError compute_hydrothermal_nusselt(10.0, 100.0, 1.1)
        @test_throws DomainError compute_hydrothermal_nusselt(
            10.0, 100.0, 0.5; phi_start=0.6, phi_end=0.4
        )
    end

    @testset "Effective Thermal Conductivity Enhancement and Damping" begin
        k_cond = 2.5

        # 1. No enhancement when Nu = 1.0
        k_eff_1 = compute_effective_hydrothermal_conductivity(k_cond, 1.0)
        @test k_eff_1 ≈ k_cond rtol = 1.0e-12

        # 2. Direct enhancement with Nu = 4.0 and picard_damping = 1.0
        k_eff_4 = compute_effective_hydrothermal_conductivity(
            k_cond, 4.0; picard_damping=1.0, resolution_weighting=false
        )
        @test k_eff_4 ≈ 10.0 rtol = 1.0e-12

        # 3. Picard relaxation damping
        k_eff_damped = compute_effective_hydrothermal_conductivity(
            k_cond, 4.0; picard_damping=0.5, k_prev=k_cond, resolution_weighting=false
        )
        @test k_eff_damped ≈ 6.25 rtol = 1.0e-12

        # 4. Cell-Péclet resolution weighting
        k_res_0 = compute_effective_hydrothermal_conductivity(
            k_cond, 4.0; Pe_cell=0.0, Pe_crit=2.0, resolution_weighting=true
        )
        @test k_res_0 ≈ 10.0 rtol = 1.0e-12

        k_res_half = compute_effective_hydrothermal_conductivity(
            k_cond, 4.0; Pe_cell=1.0, Pe_crit=2.0, resolution_weighting=true
        )
        @test k_res_half ≈ 6.25 rtol = 1.0e-12

        k_res_full = compute_effective_hydrothermal_conductivity(
            k_cond, 4.0; Pe_cell=2.5, Pe_crit=2.0, resolution_weighting=true
        )
        @test k_res_full ≈ k_cond rtol = 1.0e-12

        # 5. Clamping to floor and ceiling
        k_clamped_high = compute_effective_hydrothermal_conductivity(
            k_cond, 1.0e7; k_cutoff=500.0, resolution_weighting=false
        )
        @test k_clamped_high ≈ 500.0 rtol = 1.0e-12

        k_floor_test = compute_effective_hydrothermal_conductivity(
            k_cond, 1.0; k_floor=0.01
        )
        @test k_floor_test ≈ k_cond rtol = 1.0e-12

        # Error guards
        @test_throws DomainError compute_effective_hydrothermal_conductivity(-k_cond, 2.0)
        @test_throws DomainError compute_effective_hydrothermal_conductivity(k_cond, 0.8)
        @test_throws DomainError compute_effective_hydrothermal_conductivity(
            k_cond, 2.0; picard_damping=0.0
        )
        @test_throws DomainError compute_effective_hydrothermal_conductivity(
            k_cond, 2.0; picard_damping=1.5
        )
        @test_throws DomainError compute_effective_hydrothermal_conductivity(
            k_cond, 2.0; Pe_cell=-1.0
        )
    end

    @testset "Full Hydrothermal Convection Closure Application" begin
        cfg = HydrothermalConfig(;
            active=true,
            T_surface_ref=273.15,
            dT_min=5.0,
            picard_damping=1.0,
            resolution_weighting=false,
        )

        k_cond = 2.0
        tm_crust = 2
        tm_air = 3

        # Inactive config returns conductive baseline
        cfg_inactive = HydrothermalConfig(; active=false)
        @test apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.4, tm_crust; cfg=cfg_inactive
        ) ≈ k_cond rtol = 1.0e-12
        @test apply_hydrothermal_convection_closure(k_cond, 350.0, 0.4, tm_air; cfg=cfg) ≈
            k_cond rtol = 1.0e-12

        # Sub-freezing ice returns baseline (no fluid mobility)
        @test apply_hydrothermal_convection_closure(
            k_cond, 250.0, 0.4, tm_crust; cfg=cfg, tmfluidphase_val=273.15
        ) ≈ k_cond rtol = 1.0e-12

        # Cold surface (T <= T_surface_ref) returns baseline
        @test apply_hydrothermal_convection_closure(
            k_cond, 273.15, 0.4, tm_crust; cfg=cfg
        ) ≈ k_cond rtol = 1.0e-12

        # Active hydrothermal enhancement in warm crust aquifer
        T_warm = 350.0
        k_eff_warm = apply_hydrothermal_convection_closure(
            k_cond, T_warm, 0.35, tm_crust; cfg=cfg
        )
        @test k_eff_warm > k_cond
        @test isfinite(k_eff_warm)

        # Temperature regularization across dT_min
        k_eff_small_dT = apply_hydrothermal_convection_closure(
            k_cond, 274.15, 0.35, tm_crust; cfg=cfg
        )
        @test k_eff_small_dT >= k_cond
        @test k_eff_small_dT < k_eff_warm

        # Discriminate quadratic temperature ramp: w_T = (dT / dT_min)^2
        # dT_min = 5.0 K, T_surface_ref = 273.15 K
        # At dT = 2.5 K: w_T = (2.5 / 5.0)^2 = 0.25
        # At dT = 5.0 K: w_T = 1.0 (full Nu)
        k_at_2_5 = apply_hydrothermal_convection_closure(
            k_cond, 273.15 + 2.5, 0.35, tm_crust; cfg=cfg
        )
        k_at_5_0 = apply_hydrothermal_convection_closure(
            k_cond, 273.15 + 5.0, 0.35, tm_crust; cfg=cfg
        )
        @test k_at_2_5 >= k_cond
        @test k_at_2_5 < k_at_5_0

        # Layer thickness scaling with H_eff
        k_h5k = apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.35, tm_crust; cfg=cfg, H_eff=5000.0
        )
        k_h10k = apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.35, tm_crust; cfg=cfg, H_eff=10000.0
        )
        @test k_h10k > k_h5k > k_cond

        # Error handling - warm, cold, and inactive states
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 350.0, -0.1, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            -k_cond, 350.0, 0.3, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 200.0, -0.1, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            -k_cond, 200.0, 0.3, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 200.0, NaN, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, -10.0, 0.3, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, NaN, 0.3, tm_crust; cfg=cfg
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 200.0, -0.1, tm_crust; cfg=cfg_inactive
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.3, tm_crust; cfg=cfg, H_eff=-500.0
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.3, tm_crust; cfg=cfg, H_eff=0.0
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.3, tm_crust; cfg=cfg, H_eff=NaN
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 350.0, 0.3, tm_crust; cfg=cfg, H_eff=Inf
        )
    end

    @testset "Marker Property Integration" begin
        marknum = 4
        tm = Int32[2, 2, 2, 3]
        tkm = Float64[350.0, 260.0, 350.0, 350.0]
        phim = Float64[0.35, 0.35, 0.35, 0.99]
        XWsolidm0 = zeros(Float64, marknum)

        rhototalm = zeros(Float64, marknum)
        rhocptotalm = zeros(Float64, marknum)
        etatotalm = zeros(Float64, marknum)
        hrtotalm = zeros(Float64, marknum)
        ktotalm_base = zeros(Float64, marknum)
        ktotalm_hydro = zeros(Float64, marknum)
        tkm_rhocptotalm = zeros(Float64, marknum)
        etafluidcur_inv_kphim = zeros(Float64, marknum)
        hrsolidm = Float64[1.0e-7, 1.0e-7, 1.0e-7]
        hrfluidm = Float64[0.0, 0.0, 0.0]

        # 1. Evaluate baseline without hydrothermal convection
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm_base,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                1;
                hydrothermal_active=false,
            )
        end

        # 2. Evaluate with hydrothermal convection active
        cfg_hydro = HydrothermalConfig(; active=true, picard_damping=1.0)
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm_hydro,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                1;
                hydrothermal_active=true,
                hydrothermal_cfg=cfg_hydro,
            )
        end

        # Marker 1 (warm crust at 350 K): enhanced ktotalm
        @test ktotalm_hydro[1] > ktotalm_base[1]
        @test isfinite(ktotalm_hydro[1])

        # Marker 2 (sub-freezing crust at 260 K): identical to baseline
        @test ktotalm_hydro[2] ≈ ktotalm_base[2] rtol = 1.0e-12

        # Marker 4 (sticky air at 350 K): identical to baseline
        @test ktotalm_hydro[4] ≈ ktotalm_base[4] rtol = 1.0e-12
    end

    @testset "Dynamic Radial Gravity" begin
        R_p = 50000.0
        g_s = 0.5

        # Center singularity avoidance
        @test compute_local_radial_gravity(0.0, R_p, g_s) ≈ 0.0 atol = 1.0e-14

        # Linear interior scaling
        @test compute_local_radial_gravity(0.5 * R_p, R_p, g_s) ≈ 0.25 rtol = 1.0e-12
        @test compute_local_radial_gravity(R_p, R_p, g_s) ≈ 0.5 rtol = 1.0e-12

        # Inverse-square exterior scaling
        @test compute_local_radial_gravity(2.0 * R_p, R_p, g_s) ≈ 0.125 rtol = 1.0e-12

        # Domain error guards
        @test_throws DomainError compute_local_radial_gravity(-10.0, R_p, g_s)
        @test_throws DomainError compute_local_radial_gravity(100.0, 0.0, g_s)
        @test_throws DomainError compute_local_radial_gravity(100.0, -100.0, g_s)
        @test_throws DomainError compute_local_radial_gravity(100.0, R_p, -0.1)
    end

    @testset "Dynamic Cell Péclet Number" begin
        # Stagnant fluid gives zero Péclet number
        @test compute_cell_peclet_number(0.0, 500.0, 2.5, 1000.0, 4000.0) ≈ 0.0 atol =
            1.0e-14

        # Analytical value check:
        # kappa = 2.5 / (1000 * 4000) = 6.25e-7 m^2/s
        # Pe = 1e-6 * 500 / 6.25e-7 = 800.0
        Pe = compute_cell_peclet_number(1.0e-6, 500.0, 2.5, 1000.0, 4000.0)
        @test Pe ≈ 800.0 rtol = 1.0e-12

        # Domain error guards
        @test_throws DomainError compute_cell_peclet_number(
            -1.0, 500.0, 2.5, 1000.0, 4000.0
        )
        @test_throws DomainError compute_cell_peclet_number(
            1.0e-6, 0.0, 2.5, 1000.0, 4000.0
        )
        @test_throws DomainError compute_cell_peclet_number(
            1.0e-6, 500.0, 0.0, 1000.0, 4000.0
        )
        @test_throws DomainError compute_cell_peclet_number(1.0e-6, 500.0, 2.5, 0.0, 4000.0)
        @test_throws DomainError compute_cell_peclet_number(1.0e-6, 500.0, 2.5, 1000.0, 0.0)
    end

    @testset "Dynamic Convective Layer Thickness" begin
        R_p = 50000.0
        H_max = 10000.0
        H_min = 100.0

        # Thick shell: clamped to H_max
        H_eff_thick = compute_hydrothermal_layer_thickness(
            45000.0, R_p; R_base=10000.0, H_max=H_max, H_min=H_min
        )
        @test H_eff_thick ≈ H_max rtol = 1.0e-12

        # Thin shell: clamped to H_min
        H_eff_thin = compute_hydrothermal_layer_thickness(
            49980.0, R_p; R_base=49950.0, H_max=H_max, H_min=H_min
        )
        @test H_eff_thin ≈ H_min rtol = 1.0e-12

        # Intermediate shell: exact difference
        H_eff_mid = compute_hydrothermal_layer_thickness(
            48000.0, R_p; R_base=45000.0, H_max=H_max, H_min=H_min
        )
        @test H_eff_mid ≈ 5000.0 rtol = 1.0e-12

        # Marker inside core / impermeable basement returns H_min
        H_eff_core = compute_hydrothermal_layer_thickness(
            5000.0, R_p; R_base=10000.0, H_max=H_max, H_min=H_min
        )
        @test H_eff_core ≈ H_min rtol = 1.0e-12

        # Domain error guards
        @test_throws DomainError compute_hydrothermal_layer_thickness(
            100.0, 0.0; R_base=0.0, H_max=H_max, H_min=H_min
        )
        @test_throws DomainError compute_hydrothermal_layer_thickness(
            100.0, R_p; R_base=0.0, H_max=0.0, H_min=H_min
        )
        @test_throws DomainError compute_hydrothermal_layer_thickness(
            100.0, R_p; R_base=0.0, H_max=H_max, H_min=0.0
        )
        @test_throws DomainError compute_hydrothermal_layer_thickness(
            100.0, R_p; R_base=0.0, H_max=500.0, H_min=1000.0
        )
        @test_throws DomainError compute_hydrothermal_layer_thickness(
            100.0, R_p; R_base=60000.0, H_max=H_max, H_min=H_min
        )
    end

    @testset "Closure Coupling with Dynamic Gravity and Cell Péclet" begin
        cfg_active = HydrothermalConfig(; active=true, picard_damping=1.0)
        k_cond = 2.5
        T = 400.0
        phi = 0.35
        K_perm = 1.0e-12

        # Zero gravity: buoyancy vanishes, returns conductive conductivity exactly
        k_zero_g = apply_hydrothermal_convection_closure(
            k_cond, T, phi, 1; cfg=cfg_active, K=K_perm, gravity=0.0
        )
        @test k_zero_g ≈ k_cond rtol = 1.0e-12

        # Negative gravity rejection (both hot and sub-freezing cold markers)
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, T, phi, 1; cfg=cfg_active, K=K_perm, gravity=-0.5
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 200.0, phi, 1; cfg=cfg_active, K=K_perm, gravity=-0.5
        )

        # Unphysical convective layer thickness rejection
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, T, phi, 1; cfg=cfg_active, K=K_perm, H_eff=-500.0
        )
        @test_throws DomainError apply_hydrothermal_convection_closure(
            k_cond, 200.0, phi, 1; cfg=cfg_active, K=K_perm, H_eff=-500.0
        )

        # Non-zero gravity gives convective enhancement
        k_with_g = apply_hydrothermal_convection_closure(
            k_cond, T, phi, 1; cfg=cfg_active, K=K_perm, gravity=0.5, Pe_cell=0.0
        )
        @test k_with_g > k_cond

        # Cell Péclet damping: resolved convection suppresses subgrid enhancement
        k_damped = apply_hydrothermal_convection_closure(
            k_cond, T, phi, 1; cfg=cfg_active, K=K_perm, gravity=0.5, Pe_cell=5.0
        )
        @test k_damped ≈ k_cond rtol = 1.0e-12

        # Partial damping at intermediate Péclet number
        k_mid_pe = apply_hydrothermal_convection_closure(
            k_cond, T, phi, 1; cfg=cfg_active, K=K_perm, gravity=0.5, Pe_cell=1.0
        )
        @test k_cond < k_mid_pe < k_with_g
    end

    @testset "Marker Integration with Dynamic Hydrothermal Closure" begin
        # Set up coordinates and markers
        cfg_sim = SimulationConfig(;
            grid=GridConfig(; xsize=100000.0, ysize=100000.0),
            geometry=GeometryConfig(; rplanet=40000.0, xcenter=50000.0, ycenter=50000.0),
            hydrothermal=HydrothermalConfig(;
                active=true,
                dynamic_gravity=true,
                dynamic_layer_depth=true,
                dynamic_peclet=true,
                picard_damping=1.0,
            ),
        )
        coords = GridCoordinates(cfg_sim.grid)

        # 4 markers:
        # m = 1: at center (r = 0 m -> g = 0 -> conductive limit)
        # m = 2: in permeable crust with zero Darcy flux (Pe_cell = 0 -> full enhancement)
        # m = 3: in permeable crust with high Darcy flux (Pe_cell > Pe_crit -> subgrid damped)
        # m = 4: inside core boundary (r = 5000 m <= rcore = 10000 m -> g_hydro = 0 -> conductive limit)
        marknum = 4
        xm = [50000.0, 50000.0, 50000.0, 50000.0]
        ym = [50000.0, 75000.0, 75000.0, 55000.0]
        tm = [1, 1, 1, 1]
        tkm = [400.0, 400.0, 400.0, 400.0]
        phim = [0.35, 0.35, 0.35, 0.35]
        XWsolidm0 = [0.1, 0.1, 0.1, 0.1]
        rhototalm = zeros(marknum)
        rhocptotalm = zeros(marknum)
        etatotalm = zeros(marknum)
        hrtotalm = zeros(marknum)
        ktotalm = zeros(marknum)
        tkm_rhocptotalm = zeros(marknum)
        etafluidcur_inv_kphim = zeros(marknum)
        hrsolidm = Float64[1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7]
        hrfluidm = Float64[0.0, 0.0, 0.0, 0.0]

        # Darcy flux grids: qxD, qyD
        qxD = zeros(coords.Ny1, coords.Nx1)
        qyD = zeros(coords.Ny1, coords.Nx1)

        # Marker 3 is at x = 50000, y = 75000; set high Darcy flux in surrounding grid
        # v_Darcy ~ 1e-4 m/s -> Pe_cell >> Pe_crit
        qyD .= 1.0e-4

        # Realistic planetesimal mass: ~8.0e17 kg for R = 40 km (rho ~ 3000 kg/m^3)
        M_pl_realistic = 8.0e17

        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                1;
                hydrothermal_active=true,
                hydrothermal_cfg=cfg_sim.hydrothermal,
                xm=xm,
                ym=ym,
                coords=coords,
                rplanet_val=40000.0,
                M_planet_val=M_pl_realistic,
                rcore_val=10000.0,
                qxD_val=m == 3 ? qxD : zeros(coords.Ny1, coords.Nx1),
                qyD_val=m == 3 ? qyD : zeros(coords.Ny1, coords.Nx1),
            )
        end

        # Marker 1 (at center, g = 0): conductive limit
        @test ktotalm[1] < 3.0

        # Marker 2 (crust, stagnant fluid): enhanced conductivity
        @test ktotalm[2] > ktotalm[1]

        # Marker 3 (crust, resolved Darcy flux): subgrid enhancement damped
        @test ktotalm[3] < ktotalm[2]
        @test ktotalm[3] ≈ ktotalm[1] rtol = 0.05

        # Marker 4 (inside core, r = 5000 m <= rcore = 10000 m): suppressed to conductive baseline
        @test ktotalm[4] ≈ ktotalm[1] rtol = 0.05

        # Verification that dynamic toggles suppress coupling when false
        cfg_static = HydrothermalConfig(;
            active=true,
            dynamic_gravity=false,
            dynamic_layer_depth=false,
            dynamic_peclet=false,
            gravity=0.5,
            H_layer=10000.0,
            picard_damping=1.0,
        )
        ktotalm_static = zeros(marknum)
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm_static,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                1;
                hydrothermal_active=true,
                hydrothermal_cfg=cfg_static,
                xm=xm,
                ym=ym,
                coords=coords,
                rplanet_val=40000.0,
                M_planet_val=M_pl_realistic,
                rcore_val=0.0,
                qxD_val=m == 3 ? qxD : zeros(coords.Ny1, coords.Nx1),
                qyD_val=m == 3 ? qyD : zeros(coords.Ny1, coords.Nx1),
            )
        end
        # With dynamic_gravity=false, center marker 1 uses bulk gravity g=0.5 -> enhanced
        @test ktotalm_static[1] > ktotalm[1]
        # With dynamic_peclet=false, marker 3 is NOT damped by resolved Darcy flux
        @test ktotalm_static[3] > ktotalm[3]

        # Discriminative core suppression test with substantial local gravity
        # Marker at r = 25 km in a 40 km planetesimal (local g ~ 0.02 m/s^2)
        ktotal_with_core = zeros(1)
        compute_marker_properties!(
            1,
            [1],
            [400.0],
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ktotal_with_core,
            zeros(1),
            zeros(1),
            Float64[1.0e-7],
            Float64[0.0],
            [0.35],
            [0.1],
            1;
            hydrothermal_active=true,
            hydrothermal_cfg=cfg_sim.hydrothermal,
            xm=[50000.0],
            ym=[75000.0],
            coords=coords,
            rplanet_val=40000.0,
            M_planet_val=M_pl_realistic,
            rcore_val=30000.0, # r = 25000 <= rcore = 30000
            xcenter_val=50000.0,
            ycenter_val=50000.0,
        )
        ktotal_without_core = zeros(1)
        compute_marker_properties!(
            1,
            [1],
            [400.0],
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ktotal_without_core,
            zeros(1),
            zeros(1),
            Float64[1.0e-7],
            Float64[0.0],
            [0.35],
            [0.1],
            1;
            hydrothermal_active=true,
            hydrothermal_cfg=cfg_sim.hydrothermal,
            xm=[50000.0],
            ym=[75000.0],
            coords=coords,
            rplanet_val=40000.0,
            M_planet_val=M_pl_realistic,
            rcore_val=0.0, # No core
            xcenter_val=50000.0,
            ycenter_val=50000.0,
        )
        # Inside core: suppressed to conductive baseline; Without core: strongly enhanced
        @test ktotal_with_core[1] < 3.0
        @test ktotal_without_core[1] > 5.0
        @test ktotal_without_core[1] > 2.0 * ktotal_with_core[1]

        # M_planet_val === nothing gravity fallback verification
        ktotal_m_fallback = zeros(1)
        compute_marker_properties!(
            1,
            [1],
            [400.0],
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ktotal_m_fallback,
            zeros(1),
            zeros(1),
            Float64[1.0e-7],
            Float64[0.0],
            [0.35],
            [0.1],
            1;
            hydrothermal_active=true,
            hydrothermal_cfg=cfg_sim.hydrothermal,
            xm=[50000.0],
            ym=[75000.0],
            coords=coords,
            rplanet_val=40000.0,
            M_planet_val=nothing,
            rcore_val=10000.0,
            xcenter_val=50000.0,
            ycenter_val=50000.0,
        )
        # With M_planet_val=nothing, g_surf = cfg.gravity = 0.5 > mass-derived g_surf (~0.033)
        @test ktotal_m_fallback[1] > ktotalm[2]

        # Off-center planet center verification
        ktotalm_offcenter = zeros(1)
        compute_marker_properties!(
            1,
            [1],
            [400.0],
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ktotalm_offcenter,
            zeros(1),
            zeros(1),
            Float64[1.0e-7],
            Float64[0.0],
            [0.35],
            [0.1],
            1;
            hydrothermal_active=true,
            hydrothermal_cfg=cfg_sim.hydrothermal,
            xm=[30000.0],
            ym=[30000.0],
            coords=coords,
            rplanet_val=25000.0,
            M_planet_val=M_pl_realistic,
            rcore_val=5000.0,
            xcenter_val=30000.0,
            ycenter_val=30000.0,
        )
        # Marker at true planet center (30000, 30000) has r = 0 -> g = 0 -> conductive limit
        @test ktotalm_offcenter[1] < 3.0
        @test ktotalm_offcenter[1] ≈ ktotalm[1] rtol = 0.05

        # Staggered flux face averaging verification:
        # Non-uniform faces (2e-4, 0.0) average to 1e-4, matching uniform faces (1e-4, 1e-4)
        qxD_staggered = zeros(coords.Ny1, coords.Nx1)
        qyD_staggered = zeros(coords.Ny1, coords.Nx1)
        j_m2 = clamp(Int(floor(50000.0 / coords.dx)) + 1, 1, coords.Nx)
        i_m2 = clamp(Int(floor(75000.0 / coords.dy)) + 1, 1, coords.Ny)
        qxD_staggered[i_m2, j_m2] = 2.0e-4
        qxD_staggered[i_m2, max(1, j_m2 - 1)] = 0.0

        qxD_uniform = zeros(coords.Ny1, coords.Nx1)
        qyD_uniform = zeros(coords.Ny1, coords.Nx1)
        qxD_uniform[i_m2, j_m2] = 1.0e-4
        qxD_uniform[i_m2, max(1, j_m2 - 1)] = 1.0e-4

        ktotalm_staggered = zeros(1)
        compute_marker_properties!(
            1,
            [1],
            [400.0],
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ktotalm_staggered,
            zeros(1),
            zeros(1),
            Float64[1.0e-7],
            Float64[0.0],
            [0.35],
            [0.1],
            1;
            hydrothermal_active=true,
            hydrothermal_cfg=cfg_sim.hydrothermal,
            xm=[50000.0],
            ym=[75000.0],
            coords=coords,
            rplanet_val=40000.0,
            M_planet_val=M_pl_realistic,
            rcore_val=10000.0,
            qxD_val=qxD_staggered,
            qyD_val=qyD_staggered,
            xcenter_val=50000.0,
            ycenter_val=50000.0,
        )
        ktotalm_uniform = zeros(1)
        compute_marker_properties!(
            1,
            [1],
            [400.0],
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ktotalm_uniform,
            zeros(1),
            zeros(1),
            Float64[1.0e-7],
            Float64[0.0],
            [0.35],
            [0.1],
            1;
            hydrothermal_active=true,
            hydrothermal_cfg=cfg_sim.hydrothermal,
            xm=[50000.0],
            ym=[75000.0],
            coords=coords,
            rplanet_val=40000.0,
            M_planet_val=M_pl_realistic,
            rcore_val=10000.0,
            qxD_val=qxD_uniform,
            qyD_val=qyD_uniform,
            xcenter_val=50000.0,
            ycenter_val=50000.0,
        )
        # Staggered (2e-4 + 0)/2 equals uniform 1e-4 exactly
        @test ktotalm_staggered[1] ≈ ktotalm_uniform[1] rtol = 1.0e-6
        @test ktotalm_staggered[1] < ktotalm[2]

        # Picard damping on inactive / non-convective markers
        k_damped_inactive = apply_hydrothermal_convection_closure(
            2.0,
            400.0,
            0.35,
            1;
            cfg=HydrothermalConfig(picard_damping=0.4),
            active=false,
            k_prev=10.0,
        )
        # 0.4 * 2.0 + 0.6 * 10.0 = 0.8 + 6.0 = 6.8
        @test k_damped_inactive ≈ 6.8 rtol = 1.0e-6

        # H_eff=nothing contract accepts nothing and defaults to cfg.H_layer
        k_heff_nothing = apply_hydrothermal_convection_closure(
            2.0, 400.0, 0.35, 1; cfg=HydrothermalConfig(H_layer=12000.0), H_eff=nothing
        )
        k_heff_explicit = apply_hydrothermal_convection_closure(
            2.0, 400.0, 0.35, 1; cfg=HydrothermalConfig(H_layer=12000.0), H_eff=12000.0
        )
        @test k_heff_nothing ≈ k_heff_explicit rtol = 1.0e-6
    end
end
