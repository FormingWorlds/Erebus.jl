"""
Configuration module for Erebus.jl simulation parameters.

Provides structured types, TOML parsing, physical bounds validation,
and serialization for parameter space exploration.
"""

using DocStringExtensions
using StaticArrays
using TOML

"""
Grid discretization and spatial domain parameters.

$(FIELDS)
"""
Base.@kwdef struct GridConfig
    xsize::Float64 = 140_000.0
    ysize::Float64 = 140_000.0
    Nx::Int = 33
    Ny::Int = 33
end

"""
Planetary body geometry parameters.

$(FIELDS)
"""
Base.@kwdef struct GeometryConfig
    rplanet::Float64 = 50_000.0
    rcrust::Float64 = 50_000.0
    xcenter::Float64 = 70_000.0
    ycenter::Float64 = 70_000.0
    psurface::Float64 = 1.0e+3
    spherical_metric::Bool = false
    metric_regularization_cells::Float64 = 0.5
end

"""
Timestepping and temporal integration parameters.

Time parameters (`dt_initial`, `dt_longest`, `start_time`, `endtime`) are defined in years [yr].
Simulation routines convert them to seconds using `yearlength`.

$(FIELDS)
"""
struct TimeConfig
    dt_initial::Float64
    dt_longest::Float64
    dtcoefdn::Float64
    dtcoefup::Float64
    dtstep::Int
    dxymax::Float64
    vpratio::Float64
    DTmax::Float64
    yearlength::Float64
    start_time::Float64
    endtime::Float64
    start_step::Int
    n_steps::Int
end

function TimeConfig(;
    dt_initial::Real=1.0e11 / (365.25 * 24 * 3600),
    dt_longest::Real=1.0e11 / (365.25 * 24 * 3600),
    dtcoefdn::Real=0.5,
    dtcoefup::Real=1.2,
    dtstep::Integer=200,
    dxymax::Real=0.05,
    vpratio::Real=1.0 / 3.0,
    DTmax::Real=20.0,
    yearlength::Real=365.25 * 24 * 3600,
    start_time::Real=2.25e6,
    endtime::Real=15.0e6,
    start_step::Integer=1,
    n_steps::Integer=10,
)
    return TimeConfig(
        Float64(dt_initial),
        Float64(dt_longest),
        Float64(dtcoefdn),
        Float64(dtcoefup),
        Int(dtstep),
        Float64(dxymax),
        Float64(vpratio),
        Float64(DTmax),
        Float64(yearlength),
        Float64(start_time),
        Float64(endtime),
        Int(start_step),
        Int(n_steps),
    )
end

"""
Nonlinear iterations and solver control parameters.

$(FIELDS)
"""
Base.@kwdef struct SolverConfig
    titermax::Int = 10_000
    nplast::Int = 100_000
    yerrmax::Float64 = 1.0e+2
    etawt::Float64 = 0.0
    dphimax::Float64 = 100.01
    seed::Int = 42
    use_pardiso::Bool = false
    etaphikoef::Float64 = 1.0
    etamin::Float64 = 1.0e+12
    etamax::Float64 = 1.0e+23
    p2m_mode::Symbol = :tiled
    tile_size::Int = 4
    hydromech_solver::Symbol = :direct
    krylov_method::Symbol = :fgmres
    krylov_rtol::Float64 = 1.0e-6
    krylov_atol::Float64 = 1.0e-10
    krylov_maxiter::Int = 200
    krylov_restart::Int = 50
    darcy_elimination::Bool = false
    preconditioner::Symbol = :block_schur
    mg_levels::Int = 4
    mg_pre_smooth::Int = 2
    mg_post_smooth::Int = 2
    mg_smoother::Symbol = :damped_jacobi
    mg_omega::Float64 = 0.67
end

"""
MPI distributed execution configuration.

$(FIELDS)
"""
Base.@kwdef struct MPIConfig
    enable::Bool = false
    px::Int = 0
    py::Int = 0
    halo_width::Int = 1
end

"""
Poroelastic constitutive parameters and porosity limits.

Default compressibilities default to 0.0 to match the test baseline in `constants.jl`.
Production runs should set `betasolid = 2.5e-11` and `betafluid = 4.0e-10`.

$(FIELDS)
"""
Base.@kwdef struct PoroelasticConfig
    betasolid::Float64 = 0.0
    betafluid::Float64 = 0.0
    phimin::Float64 = 1.0e-4
    phimax::Float64 = 0.9999
    hydrofracture::Bool = false
    kappa_frac::Float64 = 1.0e3
    gamma_frac::Float64 = 1.0
    k_frac_max::Float64 = 1.0e-9
end

"""
Thermodynamic, radioactive heating, and phase change parameters.

$(FIELDS)
"""
Base.@kwdef struct ThermalConfig
    hr_al::Bool = true
    hr_fe::Bool = false
    ratio_al::Float64 = 5.0e-5
    E_al::Float64 = 5.0470e-13
    f_al::Float64 = 1.9e23
    t_half_al::Float64 = 717_000.0 * 31_540_000.0
    ratio_fe::Float64 = 1.0e-6
    E_fe::Float64 = 4.34e-13
    f_fe::Float64 = 1.957e24
    t_half_fe::Float64 = 2_620_000.0 * 31_540_000.0
    tmsolidphase::Float64 = 1416.0
    tmfluidphase::Float64 = 273.0
    Lᶠ::Float64 = 333.55e3
    phim0::Float64 = 0.2
    thermal_buoyancy::Bool = true
    fluid_viscosity_mode::Symbol = :arrhenius
    fluid_viscosity_Ea::Float64 = 15.0e3
    fluid_viscosity_T0::Float64 = 293.15
    fluid_viscosity_eta0::Float64 = 1.0e-3
    surface_radiation::Bool = false
    emissivity::Float64 = 0.9
    sigma_sb::Float64 = 5.670374419e-8
end

"""
Material phase properties for 3-phase staggered grid.

Index 1: Planetesimal core / mantle
Index 2: Porous silicate crust / rock
Index 3: Sticky air / space

$(FIELDS)
"""
Base.@kwdef struct MaterialConfig
    rhosolidm::SVector{3,Float64} = SVector{3,Float64}([3300.0, 3300.0, 1.0])
    rhofluidm::SVector{3,Float64} = SVector{3,Float64}([1000.0, 1000.0, 1.0])
    etasolidm::SVector{3,Float64} = SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
    etasolidmm::SVector{3,Float64} = SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
    etafluidm::SVector{3,Float64} = SVector{3,Float64}([1.0e+12, 1.0e+12, 1.0e-03])
    etafluidmm::SVector{3,Float64} = SVector{3,Float64}([1.0e-03, 1.0e-03, 1.0e-03])
    rhocpsolidm::SVector{3,Float64} = SVector{3,Float64}([3.3e+06, 3.3e+06, 3.0e+06])
    rhocpfluidm::SVector{3,Float64} = SVector{3,Float64}([1.0e+06, 1.0e+06, 3.0e+06])
    alphasolidm::SVector{3,Float64} = SVector{3,Float64}([3.0e-05, 3.0e-05, 0.0])
    alphafluidm::SVector{3,Float64} = SVector{3,Float64}([5.0e-05, 5.0e-05, 0.0])
    ksolidm::SVector{3,Float64} = SVector{3,Float64}([3.0, 3.0, 3000.0])
    kfluidm::SVector{3,Float64} = SVector{3,Float64}([50.0, 50.0, 3000.0])
    gggsolidm::SVector{3,Float64} = SVector{3,Float64}([1.0e+10, 1.0e+10, 1.0e+10])
    frictsolidm::SVector{3,Float64} = SVector{3,Float64}([0.6, 0.6, 0.0])
    cohessolidm::SVector{3,Float64} = SVector{3,Float64}([1.0e+08, 1.0e+08, 1.0e+08])
    tenssolidm::SVector{3,Float64} = SVector{3,Float64}([6.0e+07, 6.0e+07, 6.0e+07])
    kphim0::SVector{3,Float64} = SVector{3,Float64}([1.0e-13, 1.0e-13, 1.0e-17])
    tkm0::SVector{3,Float64} = SVector{3,Float64}([170.0, 170.0, 170.0])
    XWsolidm_init::SVector{3,Float64} = SVector{3,Float64}([0.5, 0.5, NaN])
end

"""
Output and storage parameters.

$(FIELDS)
"""
Base.@kwdef struct OutputConfig
    output_dir::String = "output"
    savematstep::Int = 10
    visstep::Int = 1
    restart_from::String = ""
    mode::Symbol = :snapshots
    telemetrystep::Int = 1
    telemetry_file::String = "telemetry.csv"
    save_final::Bool = true
end

"""
Protoplanetary disk ambient temperature evolution parameters.

Scaling exponents calibrate semi-analytical disk temperature profiles across
orbital radius and host star mass against multi-zone disk simulation models
(Drążkowska & Dullemond 2018; Lichtenberg et al. 2021; Williams et al. 2026).

$(FIELDS)
"""
Base.@kwdef struct DiskConfig
    enabled::Bool = false
    model::Symbol = :fixed
    t_ambient::Float64 = 170.0
    orbital_distance_au::Float64 = 2.5
    stellar_mass_msun::Float64 = 1.0
    t_cloud::Float64 = 30.0
    t_irr_1au::Float64 = 150.0
    t_peak_1au::Float64 = 520.0
    t_peak_time_1au_myr::Float64 = 0.12
    t_visc_0_myr::Float64 = 0.25
    gamma::Float64 = 1.4
    alpha::Float64 = 2.0
    q_irr::Float64 = 3.0 / 7.0
    q_visc::Float64 = 0.75
    p_r_t::Float64 = 0.25
    p_m_irr::Float64 = 0.25
    p_m_visc::Float64 = 0.30
    p_m_t::Float64 = 0.40
    p_m_visc_decay::Float64 = 0.30
    t_dispersal_myr::Float64 = 3.0
    dt_dispersal_myr::Float64 = 0.1
    p_amb_disk::Float64 = 10.0
    p_amb_space::Float64 = 1.0e-4
    albedo::Float64 = 0.06
    t_eq_custom::Float64 = NaN
    dispersal_active::Bool = false
end

"""
Silicate hydration and dehydration thermochemical reaction parameters.

$(FIELDS)
"""
Base.@kwdef struct ReactionConfig
    active::Bool = true
    hydration_active::Bool = true
    dehydration_active::Bool = true
    hydration_mode::Int = 1
    dehydration_mode::Int = 2
    dtreaction_hydration::Float64 = 1.0e10
    dtreaction_dehydration::Float64 = 1.0e8
    delta_H::Float64 = 40_000.0
    delta_S::Float64 = 60.0
    A_I::Float64 = 1.0e-11
    b_I::Float64 = 2.5e-4
    c_I::Float64 = 543.0
    Sxo_B::Float64 = 2.0e-11
    Tscl_B::Float64 = 10.0
    To_B::Float64 = 293.0
    alpha_relaxation::Float64 = 0.5
    pfcoeff::Float64 = 0.5
    pferrmax::Float64 = 1.0e5
    p_cavitation::Float64 = 1.0e7
    cfl_reaction::Float64 = 0.5
    dphi_reaction_max::Float64 = 0.01
    fluid_overpressure_coupling::Bool = true
end

"""
Silicate rock melting, latent heat buffering, and melt-weakened rheology parameters.

$(FIELDS)
"""
Base.@kwdef struct MeltingConfig
    active::Bool = false
    T_solidus::SVector{3,Float64} = SVector{3,Float64}([1400.0, 1400.0, NaN])
    T_liquidus::SVector{3,Float64} = SVector{3,Float64}([1800.0, 1800.0, NaN])
    L_melt::Float64 = 4.0e5
    rho_melt::Float64 = 2800.0
    alpha_eta::Float64 = 28.0
    phi_crit::Float64 = 0.4
    eta_melt::Float64 = 10.0
    dpdt_clapeyron::Float64 = 0.0
    latent_heat_mode::Symbol = :apparent_cp
    soft_turbulence::Bool = false
    turb_exponent::Float64 = 1.0 / 3.0
    eta_fluid_silicate::Float64 = 100.0
    F_turb_start::Float64 = 0.30
    F_turb_end::Float64 = 0.50
    dT_turb_min::Float64 = 10.0
    T_surface_ref::Float64 = 300.0
    k_turb_cutoff::Float64 = 1.0e6
    k_turb_floor::Float64 = 1.0e-3
end

"""
Magma ocean multi-component volatile partitioning and degassing configuration.

Configures equilibrium and dynamic degassing of H, C, N, S volatiles between silicate melt
and overlying atmosphere across magma ocean differentiation and melt crystallization.

$(FIELDS)
"""
Base.@kwdef struct MagmaOceanDegassingConfig
    active::Bool = false
    mode::Symbol = :dynamic_flux
    F_melt_threshold::Float64 = 0.40
    degas_depth_fraction::Float64 = 0.90
    crystallization_degassing::Bool = true
    redox_coupled::Bool = true
    efficiency::Float64 = 1.0

    function MagmaOceanDegassingConfig(
        active::Bool,
        mode::Symbol,
        F_melt_threshold::Real,
        degas_depth_fraction::Real,
        crystallization_degassing::Bool,
        redox_coupled::Bool,
        efficiency::Real,
    )
        (mode === :equilibrium || mode === :dynamic_flux) ||
            throw(ArgumentError("mode must be :equilibrium or :dynamic_flux, got '$mode'"))
        (0.0 <= F_melt_threshold <= 1.0) ||
            throw(DomainError(F_melt_threshold, "F_melt_threshold must be in [0, 1]"))
        (0.0 <= degas_depth_fraction <= 1.0) || throw(
            DomainError(degas_depth_fraction, "degas_depth_fraction must be in [0, 1]")
        )
        (0.0 < efficiency <= 1.0) ||
            throw(DomainError(efficiency, "efficiency must be in (0, 1]"))
        return new(
            active,
            mode,
            Float64(F_melt_threshold),
            Float64(degas_depth_fraction),
            crystallization_degassing,
            redox_coupled,
            Float64(efficiency),
        )
    end
end

"""
Two-phase buoyant silicate melt segregation and magma migration parameters.

Configures Darcy percolation of buoyant silicate melt through a compacting solid
matrix and Stokes crystal settling through magma mush/ocean regimes, including
CFL-limited subcycling, potential energy dissipation, and subsolidus freezing.

$(FIELDS)
"""
Base.@kwdef struct MagmaTransportConfig
    active::Bool = false
    k_melt_ref::Float64 = 1.0e-11
    perm_exponent::Float64 = 3.0
    phi0::Float64 = 0.10
    phi_residual::Float64 = 0.01
    phi_crit::Float64 = 0.40
    phi_pack::Float64 = 1.0
    eta_melt::Float64 = 10.0
    r_grain::Float64 = 1.0e-3
    hindered_exponent::Float64 = 2.0
    F_perc_end::Float64 = 0.35
    F_settle_start::Float64 = 0.45
    cfl_melt::Float64 = 0.5
    max_subcycles::Int = 2000
    segregation_heating::Bool = true
    latent_crystallization::Bool = true
    exsolution_active::Bool = true
    track_depletion::Bool = true
    compaction_active::Bool = false
    bulk_viscosity_ratio::Float64 = 1.0
    min_bulk_porosity::Float64 = 0.005
    compaction_length_min::Float64 = 100.0
    compaction_length_max::Float64 = 50000.0
    ponding_active::Bool = false
    eruption_active::Bool = false
    tensile_strength::Float64 = 1.0e7
    sensible_heat_transport::Bool = true
    cp_melt::Float64 = 1200.0
    sill_cooling_active::Bool = true
    crystallization_timescale::Float64 = 0.0
end

"""
Planetesimal surface volatile degassing and venting parameters.

$(FIELDS)
"""
Base.@kwdef struct VentingConfig
    active::Bool = false
    mode::Symbol = :darcy_sink
    species::Symbol = :H2O
    k_vent::Float64 = 1.0e-11
    conductance_factor::Float64 = 1.0
    L_sublimation::Float64 = 2.83e6
    latent_cooling::Bool = true
    ice_sealing::Bool = false
    t_freeze::Float64 = 273.15
    dt_seal::Float64 = 10.0
    k_seal_min_ratio::Float64 = 1.0e-6
end

"""
Multi-species volatile solubility and organic devolatilization parameters.

Configures thermodynamic speciation and multi-species volatile solubility routines
in `Erebus.Physics`. Couples marker volatile exsolution to 10-species atmospheric accumulation and kinetic escape.

$(FIELDS)
"""
Base.@kwdef struct VolatilesConfig
    active::Bool = false
    speciation_active::Bool = false
    fO2_delta_IW::Float64 = -1.0
    water_solubility_coeff::Float64 = 0.40
    water_law::Symbol = :burnham_dixon
    h2_active::Bool = false
    h2_law::Symbol = :hirschmann2012
    nitrogen_law::Symbol = :dasgupta2022
    nitrogen_henry_coeff::Float64 = 0.40
    nitrogen_nitride_capacity::Float64 = 1.0e-3
    t_organic_devol::Float64 = 550.0
    dt_organic_devol::Float64 = 50.0
    organic_n_initial_ppm::Float64 = 500.0
    initial_water_wtpct::Float64 = 1.0
    initial_carbon_ppm::Float64 = 500.0
    initial_nitrogen_ppm::Float64 = 50.0
    initial_sulfur_ppm::Float64 = 1000.0

    # Carbon solubility parameters
    carbon_active::Bool = false
    co_law::Symbol = :armstrong2015
    ch4_law::Symbol = :ardia2013
    co2_law::Symbol = :dixon1995
    graphite_saturation::Bool = true

    # Sulfur solubility parameters
    sulfur_active::Bool = false
    sulfide_law::Symbol = :boulliung2023
    sulfide_melt::Symbol = :basalt
    include_sulfate::Bool = false
    scss_active::Bool = true
    scss_law::Symbol = :smythe2017
    melt_feo_wtpct::Float64 = 10.0

    # Silicate melt composition mole fractions
    x_sio2::Float64 = 0.56
    x_al2o3::Float64 = 0.11
    x_tio2::Float64 = 0.01
end

"""
Thermodynamic volatile retention floor in nominally anhydrous minerals (NAMs) and refractory phases,
and coupling of hydrothermal surface venting to dissolved volatile depletion.

$(FIELDS)
"""
Base.@kwdef struct RetentionConfig
    active::Bool = false
    h2o_retention_ppm::Float64 = 50.0
    carbon_retention_ppm::Float64 = 50.0
    nitrogen_retention_ppm::Float64 = 5.0
    sulfur_retention_ppm::Float64 = 100.0
    T_solidus_ref::Float64 = 1400.0
    dT_retention::Float64 = 200.0
    retention_law::Symbol = :nams_exponential
    venting_drainage_active::Bool = true
    chi_vent::Float64 = 1.0
end

"""
Planetesimal redox state and electron budget configuration (Evans 2012).

Tracks extensive redox budget RB [mol e-] and dynamically buffers oxygen
fugacity relative to iron-wustite (IW) across planetary differentiation and venting.

$(FIELDS)
"""
Base.@kwdef struct RedoxConfig
    active::Bool = false
    reference::Symbol = :mantle
    serpentinization_redox::Bool = true
    segregation_redox::Bool = true
    venting_redox::Bool = true
    deltaIW_min::Float64 = -6.0
    deltaIW_max::Float64 = 6.0
    initial_x_ferric::Float64 = 0.05

    function RedoxConfig(
        active::Bool,
        reference::Symbol,
        serpentinization_redox::Bool,
        segregation_redox::Bool,
        venting_redox::Bool,
        deltaIW_min::Real,
        deltaIW_max::Real,
        initial_x_ferric::Real,
    )
        (reference === :mantle || reference === :crust) ||
            throw(ArgumentError("reference must be :mantle or :crust, got '$reference'"))
        deltaIW_min <= deltaIW_max || throw(
            ArgumentError(
                "deltaIW_min ($deltaIW_min) must be <= deltaIW_max ($deltaIW_max)"
            ),
        )
        (0.0 <= initial_x_ferric <= 1.0) ||
            throw(DomainError(initial_x_ferric, "initial_x_ferric must be in [0, 1]"))
        return new(
            active,
            reference,
            serpentinization_redox,
            segregation_redox,
            venting_redox,
            Float64(deltaIW_min),
            Float64(deltaIW_max),
            Float64(initial_x_ferric),
        )
    end
end

"""
Refractory carbon, nitrogen, sulfur, and phosphorus component configuration.

Grounds refractory element fractions and thermal breakdown thresholds across planetesimal
differentiation and accretion following:
- Carbon (f_refr_C): Bergin et al. (2026)
- Nitrogen (f_refr_N): Alexander et al. (2012)
- Sulfur (f_refr_S): Kama et al. (2019)
- Phosphorus (f_refr_P): Pasek (2008)
- Hydrogen (f_refr_H): Alexander et al. (2012), Hirschmann et al. (2006)

$(FIELDS)
"""
Base.@kwdef struct RefractoryConfig
    active::Bool = false
    f_refr_C::Float64 = 0.60
    f_refr_N::Float64 = 0.10
    f_refr_S::Float64 = 0.89
    f_refr_P::Float64 = 0.98
    f_refr_H::Float64 = 0.05
    T_pyrolysis_C::Float64 = 600.0
    T_dehydrate_H::Float64 = 750.0
    kinetics_active::Bool = false
    A_C::Float64 = 1.0e14
    Ea_C::Float64 = 2.0e5
    A_N::Float64 = 1.0e14
    Ea_N::Float64 = 2.0e5
    A_H::Float64 = 1.0e14
    Ea_H::Float64 = 1.8e5
    dh_pyro_C::Float64 = 5.0e5
    dh_pyro_N::Float64 = 5.0e5
    dh_pyro_H::Float64 = 5.0e5
    T_pyro_min::Float64 = 300.0

    function RefractoryConfig(
        active::Bool,
        f_refr_C::Real,
        f_refr_N::Real,
        f_refr_S::Real,
        f_refr_P::Real,
        f_refr_H::Real,
        T_pyrolysis_C::Real,
        T_dehydrate_H::Real,
        kinetics_active::Bool,
        A_C::Real,
        Ea_C::Real,
        A_N::Real,
        Ea_N::Real,
        A_H::Real,
        Ea_H::Real,
        dh_pyro_C::Real,
        dh_pyro_N::Real,
        dh_pyro_H::Real,
        T_pyro_min::Real,
    )
        (0.0 <= f_refr_C <= 1.0) ||
            throw(DomainError(f_refr_C, "f_refr_C must be in [0, 1]"))
        (0.0 <= f_refr_N <= 1.0) ||
            throw(DomainError(f_refr_N, "f_refr_N must be in [0, 1]"))
        (0.0 <= f_refr_S <= 1.0) ||
            throw(DomainError(f_refr_S, "f_refr_S must be in [0, 1]"))
        (0.0 <= f_refr_P <= 1.0) ||
            throw(DomainError(f_refr_P, "f_refr_P must be in [0, 1]"))
        (0.0 <= f_refr_H <= 1.0) ||
            throw(DomainError(f_refr_H, "f_refr_H must be in [0, 1]"))
        T_pyrolysis_C >= 0.0 ||
            throw(DomainError(T_pyrolysis_C, "T_pyrolysis_C must be >= 0"))
        T_dehydrate_H >= 0.0 ||
            throw(DomainError(T_dehydrate_H, "T_dehydrate_H must be >= 0"))
        T_pyro_min >= 0.0 || throw(DomainError(T_pyro_min, "T_pyro_min must be >= 0"))
        if kinetics_active
            A_C > 0.0 || throw(DomainError(A_C, "A_C must be positive"))
            Ea_C > 0.0 || throw(DomainError(Ea_C, "Ea_C must be positive"))
            A_N > 0.0 || throw(DomainError(A_N, "A_N must be positive"))
            Ea_N > 0.0 || throw(DomainError(Ea_N, "Ea_N must be positive"))
            A_H > 0.0 || throw(DomainError(A_H, "A_H must be positive"))
            Ea_H > 0.0 || throw(DomainError(Ea_H, "Ea_H must be positive"))
            dh_pyro_C >= 0.0 ||
                throw(DomainError(dh_pyro_C, "dh_pyro_C must be non-negative"))
            dh_pyro_N >= 0.0 ||
                throw(DomainError(dh_pyro_N, "dh_pyro_N must be non-negative"))
            dh_pyro_H >= 0.0 ||
                throw(DomainError(dh_pyro_H, "dh_pyro_H must be non-negative"))
        end
        return new(
            active,
            Float64(f_refr_C),
            Float64(f_refr_N),
            Float64(f_refr_S),
            Float64(f_refr_P),
            Float64(f_refr_H),
            Float64(T_pyrolysis_C),
            Float64(T_dehydrate_H),
            kinetics_active,
            Float64(A_C),
            Float64(Ea_C),
            Float64(A_N),
            Float64(Ea_N),
            Float64(A_H),
            Float64(Ea_H),
            Float64(dh_pyro_C),
            Float64(dh_pyro_N),
            Float64(dh_pyro_H),
            Float64(T_pyro_min),
        )
    end
end

"""
Multi-species volatile ice and pore fluid mixture configuration.

Parameterizes multi-snowline disk condensation and composition-dependent
freezing point depression in the H-C-N-S-P-O volatile system.

$(FIELDS)
"""
Base.@kwdef struct VolatileMixtureConfig
    active::Bool = false
    X_ice_H2O::Float64 = 0.85
    X_ice_CO2::Float64 = 0.08
    X_ice_CO::Float64 = 0.02
    X_ice_CH4::Float64 = 0.01
    X_ice_NH3::Float64 = 0.03
    X_ice_N2::Float64 = 0.005
    X_ice_H2S::Float64 = 0.005
    X_ice_PH3::Float64 = 0.0
    T_eutectic_ammonia::Float64 = 176.0
    lambda_nh3_depression::Float64 = (273.15 - 176.0) / 0.33
    lambda_solute_depression::Float64 = 50.0
    T_freeze_floor::Float64 = 176.0
    T_cond_H2O::Float64 = 160.0
    T_cond_NH3::Float64 = 135.0
    T_cond_CO2::Float64 = 75.0
    T_cond_H2S::Float64 = 75.0
    T_cond_CH4::Float64 = 45.0
    T_cond_CO::Float64 = 25.0
    T_cond_N2::Float64 = 18.0
    T_cond_PH3::Float64 = 40.0
    alpha_P::Float64 = 0.0
    P_ref::Float64 = 1.0

    function VolatileMixtureConfig(
        active::Bool,
        X_ice_H2O::Real,
        X_ice_CO2::Real,
        X_ice_CO::Real,
        X_ice_CH4::Real,
        X_ice_NH3::Real,
        X_ice_N2::Real,
        X_ice_H2S::Real,
        X_ice_PH3::Real,
        T_eutectic_ammonia::Real,
        lambda_nh3_depression::Real,
        lambda_solute_depression::Real,
        T_freeze_floor::Real,
        T_cond_H2O::Real,
        T_cond_NH3::Real,
        T_cond_CO2::Real,
        T_cond_H2S::Real,
        T_cond_CH4::Real,
        T_cond_CO::Real,
        T_cond_N2::Real,
        T_cond_PH3::Real,
        alpha_P::Real,
        P_ref::Real,
    )
        for (name, val) in [
            ("X_ice_H2O", X_ice_H2O),
            ("X_ice_CO2", X_ice_CO2),
            ("X_ice_CO", X_ice_CO),
            ("X_ice_CH4", X_ice_CH4),
            ("X_ice_NH3", X_ice_NH3),
            ("X_ice_N2", X_ice_N2),
            ("X_ice_H2S", X_ice_H2S),
            ("X_ice_PH3", X_ice_PH3),
        ]
            (0.0 <= val <= 1.0) || throw(DomainError(val, "$name must be in [0, 1]"))
        end
        T_eutectic_ammonia >= 0.0 ||
            throw(DomainError(T_eutectic_ammonia, "T_eutectic_ammonia must be >= 0"))
        T_freeze_floor >= 0.0 ||
            throw(DomainError(T_freeze_floor, "T_freeze_floor must be >= 0"))
        alpha_P >= 0.0 || throw(DomainError(alpha_P, "alpha_P must be >= 0"))
        P_ref > 0.0 || throw(DomainError(P_ref, "P_ref must be > 0"))
        for (name, val) in [
            ("T_cond_H2O", T_cond_H2O),
            ("T_cond_NH3", T_cond_NH3),
            ("T_cond_CO2", T_cond_CO2),
            ("T_cond_H2S", T_cond_H2S),
            ("T_cond_CH4", T_cond_CH4),
            ("T_cond_CO", T_cond_CO),
            ("T_cond_N2", T_cond_N2),
            ("T_cond_PH3", T_cond_PH3),
        ]
            val >= 0.0 || throw(DomainError(val, "$name must be >= 0"))
        end
        if !(
            T_cond_H2O >= T_cond_NH3 &&
            T_cond_NH3 >= T_cond_CO2 &&
            T_cond_NH3 >= T_cond_H2S &&
            T_cond_CO2 >= T_cond_CH4 &&
            T_cond_H2S >= T_cond_CH4 &&
            T_cond_CH4 >= T_cond_PH3 &&
            T_cond_PH3 >= T_cond_CO &&
            T_cond_CO >= T_cond_N2
        )
            throw(
                ArgumentError(
                    "Unphysical snowline ordering: condensation temperatures must satisfy T_cond_H2O >= T_cond_NH3 >= max(T_cond_CO2, T_cond_H2S) >= min(T_cond_CO2, T_cond_H2S) >= T_cond_CH4 >= T_cond_PH3 >= T_cond_CO >= T_cond_N2",
                ),
            )
        end
        return new(
            active,
            Float64(X_ice_H2O),
            Float64(X_ice_CO2),
            Float64(X_ice_CO),
            Float64(X_ice_CH4),
            Float64(X_ice_NH3),
            Float64(X_ice_N2),
            Float64(X_ice_H2S),
            Float64(X_ice_PH3),
            Float64(T_eutectic_ammonia),
            Float64(lambda_nh3_depression),
            Float64(lambda_solute_depression),
            Float64(T_freeze_floor),
            Float64(T_cond_H2O),
            Float64(T_cond_NH3),
            Float64(T_cond_CO2),
            Float64(T_cond_H2S),
            Float64(T_cond_CH4),
            Float64(T_cond_CO),
            Float64(T_cond_N2),
            Float64(T_cond_PH3),
            Float64(alpha_P),
            Float64(P_ref),
        )
    end
end

"""
Coupled 1D atmosphere, Guillot semi-grey radiation, disk envelope capture, and hydrodynamic escape parameters.

$(FIELDS)
"""
Base.@kwdef struct AtmosphereConfig
    active::Bool = false
    mode::Symbol = :guillot
    kappa_ir_default::Float64 = 1.0e-2
    kappa_vis_default::Float64 = 1.0e-3
    opacities::Dict{Symbol,Float64} = Dict(
        :H2O => 1.0e-2,
        :CO2 => 1.0e-3,
        :CH4 => 2.0e-3,
        :CO => 1.0e-4,
        :N2 => 1.0e-5,
        :H2 => 1.0e-5,
        :NH3 => 5.0e-3,
        :H2S => 1.0e-3,
        :SO2 => 2.0e-3,
    )
    albedo::Float64 = 0.20
    gamma_guillot::Float64 = 0.10
    T_skin_floor::Float64 = 50.0
    f_rec::Float64 = 0.10
    tau_boil::Float64 = 1.0e4 * SEC_PER_YEAR
    crossover_active::Bool = true
    b_diff_ref::Float64 = 1.0e21
end

"""
Atmospheric Jeans kinetic escape, XUV photoevaporation, and volatile mass loss parameters.

$(FIELDS)
"""
Base.@kwdef struct EscapeConfig
    active::Bool = false
    M_planet::Float64 = 1.309e18
    R_planet::Float64 = 50_000.0
    T_exobase::Float64 = 200.0
    R_exobase::Float64 = 50_000.0
    species::Symbol = :H2O
    multi_species::Bool = false
    species_list::Vector{Symbol} = [:H2O, :H2, :CO, :CO2, :CH4, :N2, :NH3, :H2S, :S2, :SO2]
    gamma::Float64 = 1.4
    hydrodynamic::Bool = true
    xuv_driven::Bool = false
    epsilon_xuv::Float64 = 0.15
    F_xuv_1au_sat::Float64 = 1.361
    t_sat_yr::Float64 = 1.0e8
    beta_xuv::Float64 = 1.23
    r_xuv_ratio::Float64 = 1.0
    tidal_correction::Bool = true

    function EscapeConfig(
        active,
        M_planet,
        R_planet,
        T_exobase,
        R_exobase,
        species,
        multi_species,
        species_list,
        gamma,
        hydrodynamic,
        xuv_driven,
        epsilon_xuv,
        F_xuv_1au_sat,
        t_sat_yr,
        beta_xuv,
        r_xuv_ratio,
        tidal_correction,
    )
        if epsilon_xuv <= 0.0 || !isfinite(epsilon_xuv)
            throw(DomainError(epsilon_xuv, "epsilon_xuv must be > 0 and finite"))
        end
        if F_xuv_1au_sat < 0.0 || !isfinite(F_xuv_1au_sat)
            throw(DomainError(F_xuv_1au_sat, "F_xuv_1au_sat must be >= 0 and finite"))
        end
        if t_sat_yr <= 0.0 || !isfinite(t_sat_yr)
            throw(DomainError(t_sat_yr, "t_sat_yr must be > 0 and finite"))
        end
        if beta_xuv < 0.0 || !isfinite(beta_xuv)
            throw(DomainError(beta_xuv, "beta_xuv must be >= 0 and finite"))
        end
        if r_xuv_ratio < 1.0 || !isfinite(r_xuv_ratio)
            throw(DomainError(r_xuv_ratio, "r_xuv_ratio must be >= 1.0 and finite"))
        end
        return new(
            active,
            Float64(M_planet),
            Float64(R_planet),
            Float64(T_exobase),
            Float64(R_exobase),
            Symbol(species),
            multi_species,
            Vector{Symbol}(species_list),
            Float64(gamma),
            hydrodynamic,
            xuv_driven,
            Float64(epsilon_xuv),
            Float64(F_xuv_1au_sat),
            Float64(t_sat_yr),
            Float64(beta_xuv),
            Float64(r_xuv_ratio),
            tidal_correction,
        )
    end
end

"""
Iron core formation by metal percolation and gravitational settling parameters.

$(FIELDS)
"""
Base.@kwdef struct CoreFormationConfig
    percolation_active::Bool = false
    settling_active::Bool = false
    sulfur_fraction::Float64 = 0.31
    metal_density_mode::Symbol = :sanloup2000
    rho_metal::Float64 = 5450.0
    rho_metal_solid::Float64 = 5700.0
    L_metal::Float64 = 2.7e5
    eta_metal::Float64 = 1.0e-2
    k_metal::Float64 = 40.0
    rhocp_metal::Float64 = 4.0e6
    Xfe_bulk::Float64 = 0.20
    phi_pack::Float64 = 0.65
    T_eutectic::Float64 = 1213.0
    dT_metal::Float64 = 50.0
    k_metal_ref::Float64 = 1.0e-9
    perm_exponent::Float64 = 3.0
    phi_crit_perc::Float64 = 0.05
    phi_residual::Float64 = 0.02
    phi0::Float64 = 0.1
    droplet_size_mode::Symbol = :capillary_mean
    droplet_diameter_fixed::Float64 = 5.0e-3
    sigma_metal_silicate::Float64 = 1.0
    We_crit::Float64 = 10.0
    hindered_exponent::Float64 = 4.5
    hadamard_rybczynski::Bool = false
    F_settle_start::Float64 = 0.40
    F_perc_end::Float64 = 0.50
    segregation_heating::Bool = true
    cfl_settling::Float64 = 0.5
    max_subcycles::Int = 2000
end

"""
Metal-silicate volatile partitioning and core segregation transport parameters.

Configures thermodynamic exchange of H, C, N, and S between molten metallic iron
and silicate melt during core formation, including donor-cell advective transport
of metal-hosted volatiles and dynamic liquid metal density coupling.

$(FIELDS)
"""
Base.@kwdef struct MetalPartitionConfig
    active::Bool = false
    model_carbon::Symbol = :grewal2019
    model_nitrogen::Symbol = :grewal2019
    model_hydrogen::Symbol = :clesi2018
    model_sulfur::Symbol = :boujibar2014
    D_H_const::Float64 = 0.5
    D_C_const::Float64 = 500.0
    D_N_const::Float64 = 20.0
    D_S_const::Float64 = 200.0
    equilibration_rate::Float64 = 1.0
    dynamic_sulfur_density::Bool = true
    D_min::Float64 = 1.0e-4
    D_max::Float64 = 1.0e5
    initial_metal_h_ppm::Float64 = 0.0
    initial_metal_c_ppm::Float64 = 0.0
    initial_metal_n_ppm::Float64 = 0.0
    initial_metal_s_ppm::Float64 = 0.0
    core_radius_fraction::Float64 = 0.5
    phi_core_threshold::Float64 = 0.40
end

"""
Normative accessory mineral tracking and meteorite diagnostic parameters.

Configures sub-eutectic stoichiometric allocation of S, P, C, and N into solid
accessory phases (troilite FeS, schreibersite (Fe,Ni)3P, cohenite (Fe,Ni)3C,
graphite C, nitrides Fe4N/CrN/TiN) and residual metallic matrix, as well as
thermal dissolution across the eutectic transition (T_eutectic ≈ 1213 K).

$(FIELDS)
"""
Base.@kwdef struct PhaseTrackingConfig
    active::Bool = false
    T_eutectic::Float64 = 1213.0
    dT_transition::Float64 = 50.0
    bulk_P_ppm::Float64 = 1000.0
    schreibersite_ni_frac::Float64 = 0.25
    cohenite_carbide_max::Float64 = 0.0667
    nitride_mode::Symbol = :roaldite
    track_regional_modes::Bool = true
    r_core_norm::Float64 = 0.5
    r_mantle_norm::Float64 = 0.85
end

"""
Hydrothermal subgrid convection parameterization configuration.

Configures porous Rayleigh-Darcy convection, free-fluid Rayleigh convection,
cubic smoothstep porosity transition blending, cell-Péclet resolution weighting,
and Picard relaxation damping.

$(FIELDS)
"""
Base.@kwdef struct HydrothermalConfig
    active::Bool = false
    phi_start::Float64 = 0.30
    phi_end::Float64 = 0.70
    Ra_m_crit::Float64 = 4.0 * pi^2
    Ra_crit::Float64 = 1100.0
    c_porous::Float64 = 1.0
    c_free::Float64 = 0.088
    H_layer::Float64 = 10000.0
    dT_min::Float64 = 5.0
    k_floor::Float64 = 1.0e-3
    k_cutoff::Float64 = 1.0e6
    picard_damping::Float64 = 0.5
    resolution_weighting::Bool = true
    Pe_crit::Float64 = 2.0
    T_surface_ref::Float64 = 273.15
    gravity::Float64 = 0.5
    cp_fluid::Float64 = 4184.0
    alpha_fluid::Float64 = 2.0e-4
    k_fluid_ref::Float64 = 0.6
    rho_fluid_ref::Float64 = 1000.0
    mu_fluid_ref::Float64 = 1.0e-3
    kphi_ref::Float64 = 1.0e-13
    sill_coupling::Bool = true
end

"""
Planetesimal accretion engine configuration.

Configures Bondi and Hill pebble accretion rates, Safronov gravitational focusing,
runaway and oligarchic growth regimes, 3D-to-2D spherical geometric mapping,
impact heating, dynamic sticky-air to rock marker conversion, volatile inheritance
from protoplanetary disk snowline evolution, and 26Al radiogenic clock inheritance.

$(FIELDS)
"""
Base.@kwdef struct AccretionConfig
    active::Bool = false
    mode::Symbol = :pebble_hill
    M_initial::Float64 = 1.0e17
    R_initial::Float64 = 20000.0
    rho_bulk::Float64 = 3000.0
    M_target::Float64 = 1.0e20
    R_target::Float64 = 50000.0
    t_start_myr::Float64 = 0.0
    t_duration_myr::Float64 = 2.0
    dM_dt_constant::Float64 = 1.5e6
    dR_dt_constant::Float64 = 5.0e-10
    tau_growth_myr::Float64 = 0.5
    h_impact::Float64 = 0.5
    v_inf::Float64 = 0.0
    cp_rock::Float64 = 1000.0
    phi_accreted::Float64 = 0.35
    Xfe_bulk_accreted::Float64 = 0.10
    snowline_coupling::Bool = true
    T_snowline_cond::Float64 = 160.0
    XWsolid_wet::Float64 = 0.40
    XWsolid_dry::Float64 = 0.0
    XH2O_wet_wtpct::Float64 = 10.0
    XH2O_dry_wtpct::Float64 = 0.1
    XC_accreted_ppm::Float64 = 1000.0
    XN_accreted_ppm::Float64 = 100.0
    XS_accreted_ppm::Float64 = 10000.0
    Sigma_peb_0::Float64 = 50.0
    p_peb::Float64 = 1.0
    stokes_number::Float64 = 0.05
    alpha_turbulence::Float64 = 1.0e-3
    c_hill::Float64 = 1.0
    c_bondi::Float64 = 1.0
    Sigma_pl_0::Float64 = 100.0
    v_disp_kms::Float64 = 0.1
    track_accretion_time::Bool = true
    stage1_mode::Symbol = :safronov
    stage2_mode::Symbol = :pebble_auto
    stage3_mode::Symbol = :safronov
    M_onset::Float64 = NaN
    f_onset::Float64 = 1.0
    M_iso::Float64 = NaN
    f_iso::Float64 = 0.5
    transition_smoothing::Bool = true
    transition_width::Float64 = 0.10
end

"""
Telescoping domain configuration for growth from small planetesimals to lunar mass.

$(FIELDS)
"""
Base.@kwdef struct TelescopingConfig
    active::Bool = false
    r_threshold_fraction::Float64 = 0.70
    max_telescope_levels::Int = 10
    target_radius::Float64 = 1_737_000.0
    buffer_markers_per_cell::Int = 4
end

"""
Top-level simulation configuration struct containing all parameter groups.

$(FIELDS)
"""
Base.@kwdef struct SimulationConfig
    grid::GridConfig = GridConfig()
    geometry::GeometryConfig = GeometryConfig()
    time::TimeConfig = TimeConfig()
    solver::SolverConfig = SolverConfig()
    poroelasticity::PoroelasticConfig = PoroelasticConfig()
    thermodynamics::ThermalConfig = ThermalConfig()
    reaction::ReactionConfig = ReactionConfig()
    materials::MaterialConfig = MaterialConfig()
    output::OutputConfig = OutputConfig()
    disk::DiskConfig = DiskConfig()
    melting::MeltingConfig = MeltingConfig()
    venting::VentingConfig = VentingConfig()
    volatiles::VolatilesConfig = VolatilesConfig()
    retention::RetentionConfig = RetentionConfig()
    escape::EscapeConfig = EscapeConfig()
    coreformation::CoreFormationConfig = CoreFormationConfig()
    metal_partition::MetalPartitionConfig = MetalPartitionConfig()
    phase_tracking::PhaseTrackingConfig = PhaseTrackingConfig()
    hydrothermal::HydrothermalConfig = HydrothermalConfig()
    accretion::AccretionConfig = AccretionConfig()
    telescoping::TelescopingConfig = TelescopingConfig()
    refractory::RefractoryConfig = RefractoryConfig()
    volatile_mixture::VolatileMixtureConfig = VolatileMixtureConfig()
    atmosphere::AtmosphereConfig = AtmosphereConfig()
    magma_transport::MagmaTransportConfig = MagmaTransportConfig()
    redox::RedoxConfig = RedoxConfig()
    magma_degassing::MagmaOceanDegassingConfig = MagmaOceanDegassingConfig()
    mpi::MPIConfig = MPIConfig()
end

"""
Returns the default simulation configuration matching baseline constants.

$(SIGNATURES)
"""
function default_config()::SimulationConfig
    return SimulationConfig()
end

function _config_field_name(ex)
    if ex isa Symbol
        return string(ex)
    elseif Meta.isexpr(ex, :.)
        parts = Symbol[]
        curr = ex
        while Meta.isexpr(curr, :.)
            push!(
                parts, curr.args[end] isa QuoteNode ? curr.args[end].value : curr.args[end]
            )
            curr = curr.args[1]
        end
        if curr isa Symbol && curr !== :cfg
            push!(parts, curr)
        end
        return join(reverse(string.(parts)), ".")
    else
        return string(ex)
    end
end

macro check_positive(expr, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    quote
        $v = $(esc(expr))
        $v > 0.0 || throw(ArgumentError(string($lbl, " must be > 0, got ", $v)))
    end
end

macro check_nonneg(expr, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    quote
        $v = $(esc(expr))
        $v >= 0.0 || throw(ArgumentError(string($lbl, " must be >= 0, got ", $v)))
    end
end

macro check_ge(expr, bound, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    b = gensym("b")
    quote
        $v = $(esc(expr))
        $b = $(esc(bound))
        $v >= $b || throw(ArgumentError(string($lbl, " must be >= ", $b, ", got ", $v)))
    end
end

macro check_positive_finite(expr, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    quote
        $v = $(esc(expr))
        ($v > 0.0 && isfinite($v)) ||
            throw(ArgumentError(string($lbl, " must be > 0 and finite, got ", $v)))
    end
end

macro check_nonneg_finite(expr, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    quote
        $v = $(esc(expr))
        ($v >= 0.0 && isfinite($v)) ||
            throw(ArgumentError(string($lbl, " must be >= 0 and finite, got ", $v)))
    end
end

macro check_finite(expr, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    quote
        $v = $(esc(expr))
        isfinite($v) || throw(ArgumentError(string($lbl, " must be finite")))
    end
end

macro check_unit_interval(expr, label=nothing)
    lbl = label === nothing ? _config_field_name(expr) : label
    v = gensym("v")
    quote
        $v = $(esc(expr))
        (0.0 <= $v <= 1.0 && isfinite($v)) ||
            throw(ArgumentError(string($lbl, " must be in [0, 1] and finite, got ", $v)))
    end
end

"""
Validates physical bounds and numerical consistency of a `SimulationConfig`.

$(SIGNATURES)

# Raises
- `ArgumentError` if any configuration parameter violates physical bounds or compiled grid constraints.
"""
function validate_config(cfg::SimulationConfig)
    # Grid checks: must be valid dimensions and domain bounds
    @check_ge cfg.grid.Nx 3
    @check_ge cfg.grid.Ny 3
    @check_positive cfg.grid.xsize "Domain xsize"
    @check_positive cfg.grid.ysize "Domain ysize"

    # Geometry checks
    @check_positive cfg.geometry.rplanet "Planet radius"
    @check_positive cfg.geometry.rcrust "Crust radius"
    cfg.geometry.rcrust <= cfg.geometry.rplanet ||
        throw(ArgumentError("Crust radius must be <= planet radius"))
    @check_positive_finite cfg.geometry.metric_regularization_cells
    min_dist_to_boundary = min(
        cfg.geometry.xcenter,
        cfg.grid.xsize - cfg.geometry.xcenter,
        cfg.geometry.ycenter,
        cfg.grid.ysize - cfg.geometry.ycenter,
    )
    (
        cfg.geometry.xcenter > 0.0 &&
        cfg.geometry.xcenter < cfg.grid.xsize &&
        cfg.geometry.ycenter > 0.0 &&
        cfg.geometry.ycenter < cfg.grid.ysize &&
        cfg.geometry.rplanet <= min_dist_to_boundary
    ) || throw(
        ArgumentError(
            "Planet of radius $(cfg.geometry.rplanet) at ($(cfg.geometry.xcenter), $(cfg.geometry.ycenter)) must fit entirely within domain [0, $(cfg.grid.xsize)] x [0, $(cfg.grid.ysize)]",
        ),
    )

    # Time checks
    @check_positive cfg.time.dt_initial "Initial dt"
    cfg.time.dt_longest >= cfg.time.dt_initial ||
        throw(ArgumentError("dt_longest must be >= dt_initial"))
    @check_ge cfg.time.n_steps 1
    @check_nonneg cfg.time.start_time
    cfg.time.endtime > cfg.time.start_time ||
        throw(ArgumentError("endtime must be > start_time"))
    @check_ge cfg.time.start_step 1
    @check_finite cfg.time.dt_initial
    @check_finite cfg.time.dt_longest

    # Poroelasticity checks
    @check_nonneg_finite cfg.poroelasticity.betasolid
    @check_nonneg_finite cfg.poroelasticity.betafluid
    0.0 < cfg.poroelasticity.phimin < cfg.poroelasticity.phimax < 1.0 || throw(
        ArgumentError(
            "Porosity bounds must satisfy 0 < phimin < phimax < 1, got phimin=$(cfg.poroelasticity.phimin), phimax=$(cfg.poroelasticity.phimax)",
        ),
    )
    @check_nonneg_finite cfg.poroelasticity.kappa_frac
    @check_positive_finite cfg.poroelasticity.gamma_frac
    @check_positive_finite cfg.poroelasticity.k_frac_max

    # Solver checks
    @check_ge cfg.solver.titermax 1
    @check_ge cfg.solver.nplast 1
    cfg.solver.titermax <= cfg.solver.nplast || throw(
        ArgumentError(
            "titermax ($(cfg.solver.titermax)) must be <= nplast ($(cfg.solver.nplast)) to prevent array bounds overflow in plastic convergence tracking",
        ),
    )
    @check_positive cfg.solver.etamin
    cfg.solver.etamax >= cfg.solver.etamin ||
        throw(ArgumentError("etamax must be >= etamin"))
    @check_positive cfg.solver.etaphikoef
    cfg.solver.p2m_mode in (:tiled, :buffered) || throw(
        ArgumentError(
            "solver.p2m_mode must be :tiled or :buffered, got :$(cfg.solver.p2m_mode)"
        ),
    )
    @check_ge cfg.solver.tile_size 2
    cfg.solver.hydromech_solver in (:direct, :iterative, :matrix_free) || throw(
        ArgumentError(
            "solver.hydromech_solver must be :direct, :iterative, or :matrix_free, got :$(cfg.solver.hydromech_solver)",
        ),
    )
    cfg.solver.krylov_method in (:fgmres, :gmres, :bicgstab) || throw(
        ArgumentError(
            "solver.krylov_method must be :fgmres, :gmres, or :bicgstab, got :$(cfg.solver.krylov_method)",
        ),
    )
    @check_positive_finite cfg.solver.krylov_rtol
    @check_nonneg_finite cfg.solver.krylov_atol
    @check_ge cfg.solver.krylov_maxiter 1
    @check_ge cfg.solver.krylov_restart 1
    cfg.solver.preconditioner in (:none, :diagonal, :block_schur, :multigrid) || throw(
        ArgumentError(
            "solver.preconditioner must be :none, :diagonal, :block_schur, or :multigrid, got :$(cfg.solver.preconditioner)",
        ),
    )
    @check_ge cfg.solver.mg_levels 1
    @check_ge cfg.solver.mg_pre_smooth 1
    @check_ge cfg.solver.mg_post_smooth 1
    cfg.solver.mg_smoother in (:damped_jacobi, :redblack_gauss_seidel) || throw(
        ArgumentError(
            "solver.mg_smoother must be :damped_jacobi or :redblack_gauss_seidel, got :$(cfg.solver.mg_smoother)",
        ),
    )
    @check_positive_finite cfg.solver.mg_omega
    cfg.solver.mg_omega <= 1.0 ||
        throw(ArgumentError("solver.mg_omega must be <= 1.0, got $(cfg.solver.mg_omega)"))
    if cfg.solver.preconditioner == :multigrid &&
        cfg.solver.mg_levels > 1 &&
        (isodd(cfg.grid.Nx) || isodd(cfg.grid.Ny))
        @warn "solver.preconditioner is :multigrid with mg_levels=$(cfg.solver.mg_levels), but grid.Nx=$(cfg.grid.Nx) or grid.Ny=$(cfg.grid.Ny) is odd; geometric multigrid coarsening requires even dimensions"
    end
    if cfg.solver.hydromech_solver == :matrix_free
        cfg.solver.darcy_elimination || throw(
            ArgumentError(
                "solver.hydromech_solver = :matrix_free requires solver.darcy_elimination = true",
            ),
        )
        !cfg.poroelasticity.hydrofracture || throw(
            ArgumentError(
                "solver.hydromech_solver = :matrix_free does not currently support poroelasticity.hydrofracture = true; use :direct or :iterative instead",
            ),
        )
        !cfg.venting.active || throw(
            ArgumentError(
                "solver.hydromech_solver = :matrix_free does not currently support venting.active = true; use :direct or :iterative instead",
            ),
        )
    end

    # MPI checks
    if cfg.mpi.enable
        cfg.mpi.px >= 0 || throw(ArgumentError("mpi.px must be >= 0, got $(cfg.mpi.px)"))
        cfg.mpi.py >= 0 || throw(ArgumentError("mpi.py must be >= 0, got $(cfg.mpi.py)"))
        cfg.mpi.halo_width >= 1 ||
            throw(ArgumentError("mpi.halo_width must be >= 1, got $(cfg.mpi.halo_width)"))
    end

    # Output checks
    @check_ge cfg.output.savematstep 1
    @check_ge cfg.output.visstep 1
    @check_ge cfg.output.telemetrystep 1
    cfg.output.mode in (:snapshots, :telemetry, :both) || throw(
        ArgumentError(
            "output.mode must be one of :snapshots, :telemetry, :both, got :$(cfg.output.mode)",
        ),
    )
    !isempty(cfg.output.telemetry_file) ||
        throw(ArgumentError("output.telemetry_file cannot be empty"))
    if !isempty(cfg.output.restart_from)
        isfile(cfg.output.restart_from) || throw(
            ArgumentError(
                "Specified restart_from checkpoint file does not exist: '$(cfg.output.restart_from)'",
            ),
        )
        endswith(cfg.output.restart_from, ".jld2") || throw(
            ArgumentError(
                "restart_from checkpoint file must have .jld2 extension, got '$(cfg.output.restart_from)'",
            ),
        )
    end

    # Thermodynamics checks
    @check_unit_interval cfg.thermodynamics.ratio_al
    @check_unit_interval cfg.thermodynamics.ratio_fe
    cfg.thermodynamics.tmfluidphase < cfg.thermodynamics.tmsolidphase || throw(
        ArgumentError(
            "tmfluidphase ($(cfg.thermodynamics.tmfluidphase)) must be < tmsolidphase ($(cfg.thermodynamics.tmsolidphase))",
        ),
    )
    @check_positive cfg.thermodynamics.Lᶠ
    @check_positive_finite cfg.thermodynamics.E_al
    @check_positive_finite cfg.thermodynamics.f_al
    @check_positive_finite cfg.thermodynamics.t_half_al
    @check_positive_finite cfg.thermodynamics.E_fe
    @check_positive_finite cfg.thermodynamics.f_fe
    @check_positive_finite cfg.thermodynamics.t_half_fe
    cfg.thermodynamics.fluid_viscosity_mode in Set([:arrhenius, :constant]) || throw(
        ArgumentError(
            "fluid_viscosity_mode must be :arrhenius or :constant, got $(cfg.thermodynamics.fluid_viscosity_mode)",
        ),
    )
    @check_nonneg_finite cfg.thermodynamics.fluid_viscosity_Ea
    @check_positive_finite cfg.thermodynamics.fluid_viscosity_T0
    @check_positive_finite cfg.thermodynamics.fluid_viscosity_eta0
    @check_unit_interval cfg.thermodynamics.emissivity
    @check_positive_finite cfg.thermodynamics.sigma_sb

    # Disk checks
    cfg.disk.model in Set([:fixed, :monotonic, :class1_to_class2, :class0_to_class2]) ||
        throw(
            ArgumentError(
                "disk model must be :fixed, :monotonic, :class1_to_class2, or :class0_to_class2, got $(cfg.disk.model)",
            ),
        )
    @check_positive_finite cfg.disk.t_ambient
    @check_positive_finite cfg.disk.orbital_distance_au
    @check_positive_finite cfg.disk.stellar_mass_msun
    @check_positive_finite cfg.disk.t_cloud
    @check_positive_finite cfg.disk.t_irr_1au
    @check_positive_finite cfg.disk.t_peak_1au
    @check_positive_finite cfg.disk.t_peak_time_1au_myr
    @check_positive_finite cfg.disk.t_visc_0_myr
    @check_positive_finite cfg.disk.gamma
    @check_positive_finite cfg.disk.alpha
    @check_positive_finite cfg.disk.q_irr
    @check_positive_finite cfg.disk.q_visc
    @check_nonneg_finite cfg.disk.p_r_t
    @check_nonneg_finite cfg.disk.p_m_irr
    @check_nonneg_finite cfg.disk.p_m_visc
    @check_nonneg_finite cfg.disk.p_m_t
    @check_nonneg_finite cfg.disk.p_m_visc_decay
    @check_positive_finite cfg.disk.t_dispersal_myr
    @check_positive_finite cfg.disk.dt_dispersal_myr
    @check_positive_finite cfg.disk.p_amb_disk
    @check_nonneg_finite cfg.disk.p_amb_space
    (0.0 <= cfg.disk.albedo < 1.0) && isfinite(cfg.disk.albedo) ||
        throw(ArgumentError("albedo must be in [0.0, 1.0), got $(cfg.disk.albedo)"))
    isnan(cfg.disk.t_eq_custom) ||
        (cfg.disk.t_eq_custom > 0.0 && isfinite(cfg.disk.t_eq_custom)) ||
        throw(
            ArgumentError(
                "t_eq_custom must be > 0 and finite when specified, got $(cfg.disk.t_eq_custom)",
            ),
        )

    # Reaction checks
    cfg.reaction.hydration_mode in Set([1, 2, 3, 9]) || throw(
        ArgumentError(
            "hydration_mode must be 1, 2, 3, or 9, got $(cfg.reaction.hydration_mode)"
        ),
    )
    cfg.reaction.dehydration_mode in Set([1, 2, 3, 9]) || throw(
        ArgumentError(
            "dehydration_mode must be 1, 2, 3, or 9, got $(cfg.reaction.dehydration_mode)",
        ),
    )
    @check_positive_finite cfg.reaction.dtreaction_hydration
    @check_positive_finite cfg.reaction.dtreaction_dehydration
    @check_positive_finite cfg.reaction.delta_H
    @check_positive_finite cfg.reaction.delta_S
    @check_positive_finite cfg.reaction.A_I
    @check_positive_finite cfg.reaction.b_I
    @check_positive_finite cfg.reaction.c_I
    @check_positive_finite cfg.reaction.Sxo_B
    @check_positive_finite cfg.reaction.Tscl_B
    @check_positive_finite cfg.reaction.To_B
    0.0 < cfg.reaction.alpha_relaxation <= 1.0 || throw(
        ArgumentError(
            "alpha_relaxation must be in (0, 1], got $(cfg.reaction.alpha_relaxation)"
        ),
    )
    @check_unit_interval cfg.reaction.pfcoeff
    @check_positive_finite cfg.reaction.pferrmax
    @check_nonneg_finite cfg.reaction.p_cavitation
    @check_positive_finite cfg.reaction.cfl_reaction
    0.0 < cfg.reaction.dphi_reaction_max <= 1.0 &&
    isfinite(cfg.reaction.dphi_reaction_max) ||
        throw(ArgumentError("dphi_reaction_max must be in (0, 1] and finite"))

    # Materials checks: all 18 property arrays must be positive/non-negative and finite
    for (arr, name, strictly_pos) in [
        (cfg.materials.rhosolidm, "rhosolidm", true),
        (cfg.materials.rhofluidm, "rhofluidm", true),
        (cfg.materials.etasolidm, "etasolidm", true),
        (cfg.materials.etasolidmm, "etasolidmm", true),
        (cfg.materials.etafluidm, "etafluidm", true),
        (cfg.materials.etafluidmm, "etafluidmm", true),
        (cfg.materials.rhocpsolidm, "rhocpsolidm", true),
        (cfg.materials.rhocpfluidm, "rhocpfluidm", true),
        (cfg.materials.alphasolidm, "alphasolidm", false),
        (cfg.materials.alphafluidm, "alphafluidm", false),
        (cfg.materials.ksolidm, "ksolidm", true),
        (cfg.materials.kfluidm, "kfluidm", true),
        (cfg.materials.gggsolidm, "gggsolidm", true),
        (cfg.materials.frictsolidm, "frictsolidm", false),
        (cfg.materials.cohessolidm, "cohessolidm", true),
        (cfg.materials.tenssolidm, "tenssolidm", true),
        (cfg.materials.kphim0, "kphim0", true),
        (cfg.materials.tkm0, "tkm0", true),
    ]
        all(isfinite, arr) || throw(ArgumentError("$name elements must be finite"))
        if strictly_pos
            all(x -> x > 0.0, arr) || throw(ArgumentError("$name elements must be > 0"))
        else
            all(x -> x >= 0.0, arr) || throw(ArgumentError("$name elements must be >= 0"))
        end
    end

    if cfg.materials.rhosolidm != rhosolidm ||
        cfg.materials.rhofluidm != rhofluidm ||
        cfg.materials.etasolidm != etasolidm ||
        cfg.materials.etasolidmm != etasolidmm ||
        cfg.materials.etafluidm != etafluidm ||
        cfg.materials.etafluidmm != etafluidmm ||
        cfg.materials.ksolidm != ksolidm ||
        cfg.materials.kfluidm != kfluidm
        throw(
            ArgumentError(
                "Overriding [materials] arrays at runtime is not supported because marker property assignment uses compiled constants. Recompilation is required to modify material properties.",
            ),
        )
    end

    # Melting checks
    if cfg.melting.soft_turbulence && !cfg.melting.active
        throw(
            ArgumentError(
                "Melting soft_turbulence cannot be enabled when melting active is false"
            ),
        )
    end

    if cfg.melting.active
        for idx in 1:2
            T_s = cfg.melting.T_solidus[idx]
            T_l = cfg.melting.T_liquidus[idx]
            (T_s > 0.0 && isfinite(T_s)) || throw(
                ArgumentError("Melting T_solidus[$idx] must be > 0 and finite, got $T_s"),
            )
            (T_l > 0.0 && isfinite(T_l)) || throw(
                ArgumentError("Melting T_liquidus[$idx] must be > 0 and finite, got $T_l"),
            )
            T_s < T_l || throw(
                ArgumentError(
                    "Melting T_solidus[$idx] ($T_s) must be < T_liquidus[$idx] ($T_l)"
                ),
            )
        end
        @check_positive_finite cfg.melting.L_melt
        @check_positive_finite cfg.melting.rho_melt
        @check_nonneg_finite cfg.melting.alpha_eta
        (0.0 < cfg.melting.phi_crit < 1.0 && isfinite(cfg.melting.phi_crit)) || throw(
            ArgumentError(
                "Melting phi_crit must be in (0, 1), got $(cfg.melting.phi_crit)"
            ),
        )
        @check_positive_finite cfg.melting.eta_melt
        @check_nonneg_finite cfg.melting.dpdt_clapeyron
        cfg.melting.latent_heat_mode == :apparent_cp || throw(
            ArgumentError(
                "Melting latent_heat_mode must be :apparent_cp, got $(cfg.melting.latent_heat_mode)",
            ),
        )

        if cfg.melting.soft_turbulence
            @check_positive_finite cfg.melting.turb_exponent
            @check_positive_finite cfg.melting.eta_fluid_silicate
            (0.0 <= cfg.melting.F_turb_start < cfg.melting.F_turb_end <= 1.0) || throw(
                ArgumentError(
                    "Melting F_turb bounds must satisfy 0 <= F_turb_start < F_turb_end <= 1, got [$(cfg.melting.F_turb_start), $(cfg.melting.F_turb_end)]",
                ),
            )
            @check_positive_finite cfg.melting.dT_turb_min
            @check_positive_finite cfg.melting.T_surface_ref
            (
                0.0 < cfg.melting.k_turb_floor < cfg.melting.k_turb_cutoff &&
                isfinite(cfg.melting.k_turb_cutoff)
            ) || throw(
                ArgumentError(
                    "Melting k_turb bounds must satisfy 0 < k_turb_floor < k_turb_cutoff, got floor=$(cfg.melting.k_turb_floor), cutoff=$(cfg.melting.k_turb_cutoff)",
                ),
            )
        end
    end

    # Magma transport checks
    if cfg.magma_transport.active
        cfg.melting.active || throw(
            ArgumentError("magma_transport cannot be active when melting is inactive")
        )
        mt = cfg.magma_transport
        @check_positive_finite mt.k_melt_ref
        @check_positive_finite mt.perm_exponent
        (0.0 < mt.phi0 <= 1.0 && isfinite(mt.phi0)) ||
            throw(ArgumentError("magma_transport phi0 must be in (0, 1], got $(mt.phi0)"))
        (0.0 <= mt.phi_residual < mt.phi_crit && isfinite(mt.phi_residual)) || throw(
            ArgumentError(
                "magma_transport phi_residual must satisfy 0 <= phi_residual < phi_crit, got $(mt.phi_residual)",
            ),
        )
        (0.0 < mt.phi_crit < 1.0 && isfinite(mt.phi_crit)) || throw(
            ArgumentError("magma_transport phi_crit must be in (0, 1), got $(mt.phi_crit)"),
        )
        (mt.phi_crit <= mt.phi_pack <= 1.0 && isfinite(mt.phi_pack)) || throw(
            ArgumentError(
                "magma_transport phi_pack must satisfy phi_crit <= phi_pack <= 1, got $(mt.phi_pack)",
            ),
        )
        @check_positive_finite mt.eta_melt
        @check_positive_finite mt.r_grain
        @check_nonneg_finite mt.hindered_exponent
        (
            0.0 <=
            mt.phi_residual <
            mt.F_perc_end <=
            mt.phi_crit <=
            mt.F_settle_start <=
            mt.phi_pack <=
            1.0
        ) || throw(
            ArgumentError(
                "magma_transport regime ordering must satisfy 0 <= phi_residual < F_perc_end <= phi_crit <= F_settle_start <= phi_pack <= 1, got [$(mt.phi_residual), $(mt.F_perc_end), $(mt.phi_crit), $(mt.F_settle_start), $(mt.phi_pack)]",
            ),
        )
        (0.0 < mt.cfl_melt <= 1.0 && isfinite(mt.cfl_melt)) || throw(
            ArgumentError("magma_transport cfl_melt must be in (0, 1], got $(mt.cfl_melt)"),
        )
        mt.max_subcycles >= 1 || throw(
            ArgumentError(
                "magma_transport max_subcycles must be >= 1, got $(mt.max_subcycles)"
            ),
        )
        @check_positive_finite mt.bulk_viscosity_ratio
        (0.0 < mt.min_bulk_porosity < mt.phi_crit && isfinite(mt.min_bulk_porosity)) ||
            throw(
                ArgumentError(
                    "magma_transport min_bulk_porosity must satisfy 0 < min_bulk_porosity < phi_crit, got $(mt.min_bulk_porosity)",
                ),
            )
        (
            0.0 < mt.compaction_length_min <= mt.compaction_length_max &&
            isfinite(mt.compaction_length_min) &&
            isfinite(mt.compaction_length_max)
        ) || throw(
            ArgumentError(
                "magma_transport compaction length bounds must satisfy 0 < compaction_length_min <= compaction_length_max, got min=$(mt.compaction_length_min), max=$(mt.compaction_length_max)",
            ),
        )
        @check_positive_finite mt.tensile_strength
        @check_positive_finite mt.cp_melt
        @check_nonneg_finite mt.crystallization_timescale
        if mt.eruption_active && !mt.compaction_active
            throw(
                ArgumentError(
                    "magma_transport eruption_active=true requires compaction_active=true"
                ),
            )
        end
    end

    # Venting checks
    cfg.venting.mode in Set([:darcy_sink, :hydrofracture_gated]) || throw(
        ArgumentError(
            "venting mode must be :darcy_sink or :hydrofracture_gated, got $(cfg.venting.mode)",
        ),
    )
    @check_positive_finite cfg.venting.k_vent
    @check_positive_finite cfg.venting.conductance_factor
    @check_positive_finite cfg.venting.L_sublimation
    @check_positive_finite cfg.venting.t_freeze
    @check_positive_finite cfg.venting.dt_seal
    (0.0 < cfg.venting.k_seal_min_ratio <= 1.0 && isfinite(cfg.venting.k_seal_min_ratio)) ||
        throw(
            ArgumentError(
                "k_seal_min_ratio must be in (0, 1], got $(cfg.venting.k_seal_min_ratio)"
            ),
        )
    cfg.venting.species in Set([:H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2]) ||
        throw(
            ArgumentError(
                "venting species must be one of standard volatile species, got $(cfg.venting.species)",
            ),
        )

    # Volatiles checks
    (isfinite(cfg.volatiles.fO2_delta_IW) && abs(cfg.volatiles.fO2_delta_IW) <= 50.0) ||
        throw(
            ArgumentError(
                "fO2_delta_IW must be finite and within [-50, 50], got $(cfg.volatiles.fO2_delta_IW)",
            ),
        )
    @check_positive_finite cfg.volatiles.water_solubility_coeff
    @check_positive_finite cfg.volatiles.nitrogen_henry_coeff
    @check_positive_finite cfg.volatiles.nitrogen_nitride_capacity
    @check_positive_finite cfg.volatiles.t_organic_devol
    @check_positive_finite cfg.volatiles.dt_organic_devol
    @check_nonneg_finite cfg.volatiles.organic_n_initial_ppm
    @check_nonneg_finite cfg.volatiles.initial_water_wtpct
    @check_nonneg_finite cfg.volatiles.initial_carbon_ppm
    @check_nonneg_finite cfg.volatiles.initial_nitrogen_ppm
    @check_nonneg_finite cfg.volatiles.initial_sulfur_ppm
    cfg.volatiles.water_law in
    Set([:burnham_dixon, :sossi_peridotite, :basalt_dixon, :newcombe_lunar]) || throw(
        ArgumentError(
            "water_law must be :burnham_dixon, :sossi_peridotite, :basalt_dixon, or :newcombe_lunar, got $(cfg.volatiles.water_law)",
        ),
    )
    cfg.volatiles.h2_law in Set([:hirschmann2012, :gaillard2003]) || throw(
        ArgumentError(
            "h2_law must be :hirschmann2012 or :gaillard2003, got $(cfg.volatiles.h2_law)",
        ),
    )
    cfg.volatiles.nitrogen_law in Set([:dasgupta2022, :libourel2003]) || throw(
        ArgumentError(
            "nitrogen_law must be :dasgupta2022 or :libourel2003, got $(cfg.volatiles.nitrogen_law)",
        ),
    )
    cfg.volatiles.co_law in Set([:armstrong2015, :yoshioka2019_morb]) || throw(
        ArgumentError(
            "co_law must be :armstrong2015 or :yoshioka2019_morb, got $(cfg.volatiles.co_law)",
        ),
    )
    cfg.volatiles.ch4_law in Set([:ardia2013]) ||
        throw(ArgumentError("ch4_law must be :ardia2013, got $(cfg.volatiles.ch4_law)"))
    cfg.volatiles.co2_law in Set([:dixon1995]) ||
        throw(ArgumentError("co2_law must be :dixon1995, got $(cfg.volatiles.co2_law)"))
    cfg.volatiles.sulfide_law in Set([:boulliung2023, :gaillard2022]) || throw(
        ArgumentError(
            "sulfide_law must be :boulliung2023 or :gaillard2022, got $(cfg.volatiles.sulfide_law)",
        ),
    )
    cfg.volatiles.sulfide_melt in Set([:basalt, :andesite, :trachybasalt]) || throw(
        ArgumentError(
            "sulfide_melt must be :basalt, :andesite, or :trachybasalt, got $(cfg.volatiles.sulfide_melt)",
        ),
    )
    cfg.volatiles.scss_law in Set([:smythe2017, :oneill2002]) || throw(
        ArgumentError(
            "scss_law must be :smythe2017 or :oneill2002, got $(cfg.volatiles.scss_law)"
        ),
    )
    @check_nonneg_finite cfg.volatiles.melt_feo_wtpct
    @check_unit_interval cfg.volatiles.x_sio2
    @check_unit_interval cfg.volatiles.x_al2o3
    @check_unit_interval cfg.volatiles.x_tio2

    # Escape checks
    @check_positive_finite cfg.escape.M_planet
    @check_positive_finite cfg.escape.R_planet
    if cfg.escape.active
        isapprox(cfg.escape.R_planet, cfg.geometry.rplanet; rtol=0.01) || throw(
            ArgumentError(
                "escape.R_planet ($(cfg.escape.R_planet)) must match simulated planet radius geometry.rplanet ($(cfg.geometry.rplanet)) when escape.active=true",
            ),
        )
    end
    @check_positive_finite cfg.escape.T_exobase
    (cfg.escape.R_exobase >= cfg.escape.R_planet && isfinite(cfg.escape.R_exobase)) ||
        throw(
            ArgumentError(
                "R_exobase must be >= R_planet ($(cfg.escape.R_planet)) and finite, got $(cfg.escape.R_exobase)",
            ),
        )
    @check_positive_finite cfg.escape.gamma
    cfg.escape.species in Set([:H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2]) ||
        throw(
            ArgumentError(
                "escape species must be one of :H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2, got $(cfg.escape.species)",
            ),
        )
    for sp in cfg.escape.species_list
        sp in Set([:H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2]) || throw(
            ArgumentError(
                "escape species_list elements must be one of standard volatile species, got $sp",
            ),
        )
    end
    if cfg.venting.active && cfg.escape.active
        if !cfg.escape.multi_species
            cfg.venting.species == cfg.escape.species || throw(
                ArgumentError(
                    "venting species $(cfg.venting.species) must match escape species $(cfg.escape.species) when escape.multi_species=false",
                ),
            )
        else
            cfg.venting.species in cfg.escape.species_list || throw(
                ArgumentError(
                    "venting species $(cfg.venting.species) must be in escape species_list ($(cfg.escape.species_list)) when both venting and escape are active",
                ),
            )
        end
    end
    if cfg.venting.active && cfg.volatiles.speciation_active
        cfg.venting.species == :H2O || throw(
            ArgumentError(
                "venting.species must be :H2O when volatiles.speciation_active=true (received :$(cfg.venting.species))",
            ),
        )
    end
    if cfg.escape.active && cfg.volatiles.speciation_active
        cfg.escape.multi_species || throw(
            ArgumentError(
                "escape.multi_species must be true when volatiles.speciation_active=true",
            ),
        )
        issubset(SPECIATION_SPECIES, Set(cfg.escape.species_list)) || throw(
            ArgumentError(
                "escape.species_list must contain all thermodynamic speciation species when volatiles.speciation_active=true",
            ),
        )
    end
    if cfg.escape.active &&
        cfg.escape.multi_species &&
        cfg.retention.active &&
        cfg.retention.venting_drainage_active &&
        !cfg.volatiles.speciation_active
        required_drain_species = Set([:H2O, :CO2, :N2, :H2S, cfg.venting.species])
        issubset(required_drain_species, Set(cfg.escape.species_list)) || throw(
            ArgumentError(
                "escape.species_list must contain :H2O, :CO2, :N2, :H2S, and $(cfg.venting.species) when retention.venting_drainage_active=true",
            ),
        )
    end

    # Retention checks
    if cfg.retention.active && !cfg.volatiles.active
        throw(
            ArgumentError(
                "RetentionConfig active=true requires VolatilesConfig active=true (enable [volatiles] active = true).",
            ),
        )
    end
    @check_nonneg_finite cfg.retention.h2o_retention_ppm
    @check_nonneg_finite cfg.retention.carbon_retention_ppm
    @check_nonneg_finite cfg.retention.nitrogen_retention_ppm
    @check_nonneg_finite cfg.retention.sulfur_retention_ppm
    @check_positive_finite cfg.retention.T_solidus_ref
    @check_positive_finite cfg.retention.dT_retention
    cfg.retention.retention_law in
    Set([:nams_exponential, :constant_floor, :linear_melt_blend]) || throw(
        ArgumentError(
            "retention_law must be :nams_exponential, :constant_floor, or :linear_melt_blend, got $(cfg.retention.retention_law)",
        ),
    )
    @check_unit_interval cfg.retention.chi_vent

    # Core formation validation
    if cfg.coreformation.percolation_active || cfg.coreformation.settling_active
        cf = cfg.coreformation
        (0.0 <= cf.sulfur_fraction <= 0.40 && isfinite(cf.sulfur_fraction)) || throw(
            ArgumentError(
                "sulfur_fraction must be in [0, 0.40] and finite, got $(cf.sulfur_fraction)",
            ),
        )
        cf.metal_density_mode in Set([:sanloup2000, :morard2014, :constant]) || throw(
            ArgumentError(
                "metal_density_mode must be :sanloup2000, :morard2014, or :constant, got $(cf.metal_density_mode)",
            ),
        )
        (isfinite(cf.rho_metal) && cf.rho_metal > cfg.materials.rhosolidm[1]) || throw(
            ArgumentError(
                "rho_metal ($(cf.rho_metal)) must be finite and exceed silicate rock density ($(cfg.materials.rhosolidm[1]))",
            ),
        )
        (isfinite(cf.rho_metal_solid) && cf.rho_metal_solid > cf.rho_metal) || throw(
            ArgumentError(
                "rho_metal_solid ($(cf.rho_metal_solid)) must be finite and > rho_metal ($(cf.rho_metal))",
            ),
        )
        @check_nonneg_finite cf.L_metal
        @check_positive_finite cf.eta_metal
        @check_positive_finite cf.k_metal
        @check_positive_finite cf.rhocp_metal
        (0.0 <= cf.Xfe_bulk <= cf.phi_pack) || throw(
            ArgumentError(
                "Xfe_bulk must be in [0, phi_pack] ($([0, cf.phi_pack])), got $(cf.Xfe_bulk)",
            ),
        )
        (0.0 < cf.phi_pack <= 1.0) ||
            throw(ArgumentError("phi_pack must be in (0, 1], got $(cf.phi_pack)"))
        (0.0 <= cf.phi_residual <= cf.phi_crit_perc < cf.phi_pack) || throw(
            ArgumentError(
                "phi_residual ($(cf.phi_residual)) must be <= phi_crit_perc ($(cf.phi_crit_perc)) < phi_pack ($(cf.phi_pack))",
            ),
        )
        (0.0 < cf.phi0 < 1.0) ||
            throw(ArgumentError("phi0 must be in (0, 1), got $(cf.phi0)"))
        (0.0 <= cf.F_settle_start <= cf.F_perc_end <= 1.0) || throw(
            ArgumentError(
                "F_settle_start ($(cf.F_settle_start)) must be <= F_perc_end ($(cf.F_perc_end)) in [0, 1]",
            ),
        )
        (0.0 < cf.cfl_settling <= 1.0) ||
            throw(ArgumentError("cfl_settling must be in (0, 1], got $(cf.cfl_settling)"))
        cf.max_subcycles > 0 ||
            throw(ArgumentError("max_subcycles must be > 0, got $(cf.max_subcycles)"))
        cf.droplet_size_mode in
        Set([:fixed, :capillary_mean, :bond_mean, :weber_mean, :weber_turbulent]) || throw(
            ArgumentError(
                "droplet_size_mode must be one of :fixed, :capillary_mean, :bond_mean, :weber_mean, :weber_turbulent, got $(cf.droplet_size_mode)",
            ),
        )
        @check_positive_finite cf.droplet_diameter_fixed
        @check_positive_finite cf.sigma_metal_silicate
        @check_positive_finite cf.We_crit
        @check_nonneg_finite cf.hindered_exponent
    end

    if cfg.coreformation.percolation_active
        cf = cfg.coreformation
        cf.T_eutectic < cfg.melting.T_solidus[1] || throw(
            ArgumentError(
                "T_eutectic ($(cf.T_eutectic)) must be below silicate solidus ($(cfg.melting.T_solidus[1])) for percolation",
            ),
        )
        @check_positive_finite cf.dT_metal
        @check_positive_finite cf.k_metal_ref
        @check_nonneg_finite cf.perm_exponent
    end

    if cfg.coreformation.settling_active
        cfg.melting.active || throw(
            ArgumentError(
                "settling_active requires melting.active = true for silicate melt fraction",
            ),
        )
    end

    if cfg.volatiles.active && !cfg.melting.active
        @warn "VolatilesConfig active=true without melting.active=true: silicate melt volatile exsolution occurs only when melting.active=true."
    end

    # Metal partition validation
    if cfg.metal_partition.active
        if !cfg.volatiles.active
            throw(
                ArgumentError(
                    "MetalPartitionConfig active=true requires VolatilesConfig active=true (enable [volatiles] active = true).",
                ),
            )
        end
        if !cfg.coreformation.percolation_active && !cfg.coreformation.settling_active
            @warn "MetalPartitionConfig active=true without coreformation percolation_active or settling_active: metal volatile segregation transport will remain inactive."
        end
    end
    @check_positive_finite cfg.metal_partition.D_min
    (
        cfg.metal_partition.D_max >= cfg.metal_partition.D_min &&
        isfinite(cfg.metal_partition.D_max)
    ) || throw(
        ArgumentError(
            "D_max must be >= D_min and finite, got $(cfg.metal_partition.D_max)"
        ),
    )
    @check_unit_interval cfg.metal_partition.equilibration_rate
    @check_unit_interval cfg.metal_partition.core_radius_fraction
    @check_unit_interval cfg.metal_partition.phi_core_threshold
    @check_nonneg_finite cfg.metal_partition.D_H_const
    @check_nonneg_finite cfg.metal_partition.D_C_const
    @check_nonneg_finite cfg.metal_partition.D_N_const
    @check_nonneg_finite cfg.metal_partition.D_S_const
    @check_nonneg_finite cfg.metal_partition.initial_metal_h_ppm
    @check_nonneg_finite cfg.metal_partition.initial_metal_c_ppm
    @check_nonneg_finite cfg.metal_partition.initial_metal_n_ppm
    @check_nonneg_finite cfg.metal_partition.initial_metal_s_ppm
    cfg.metal_partition.model_carbon in Set([:constant, :grewal2019, :fischer2020]) ||
        throw(
            ArgumentError(
                "model_carbon must be :constant, :grewal2019, or :fischer2020, got :$(cfg.metal_partition.model_carbon)",
            ),
        )
    cfg.metal_partition.model_nitrogen in Set([:constant, :grewal2019]) || throw(
        ArgumentError(
            "model_nitrogen must be :constant or :grewal2019, got :$(cfg.metal_partition.model_nitrogen)",
        ),
    )
    cfg.metal_partition.model_hydrogen in Set([:constant, :clesi2018]) || throw(
        ArgumentError(
            "model_hydrogen must be :constant or :clesi2018, got :$(cfg.metal_partition.model_hydrogen)",
        ),
    )
    cfg.metal_partition.model_sulfur in Set([:constant, :boujibar2014]) || throw(
        ArgumentError(
            "model_sulfur must be :constant or :boujibar2014, got :$(cfg.metal_partition.model_sulfur)",
        ),
    )

    # Phase tracking validation
    if cfg.phase_tracking.active
        (
            cfg.coreformation.percolation_active ||
            cfg.coreformation.settling_active ||
            cfg.metal_partition.active
        ) || throw(
            ArgumentError(
                "phase_tracking requires core formation (percolation_active or settling_active) or metal_partition to be active",
            ),
        )
        isapprox(cfg.phase_tracking.T_eutectic, cfg.coreformation.T_eutectic; atol=1e-3) ||
            throw(
                ArgumentError(
                    "phase_tracking.T_eutectic ($(cfg.phase_tracking.T_eutectic)) must match coreformation.T_eutectic ($(cfg.coreformation.T_eutectic))",
                ),
            )
        @check_positive_finite cfg.phase_tracking.T_eutectic
        @check_positive_finite cfg.phase_tracking.dT_transition
        @check_nonneg_finite cfg.phase_tracking.bulk_P_ppm
        @check_unit_interval cfg.phase_tracking.schreibersite_ni_frac
        (
            0.0 < cfg.phase_tracking.cohenite_carbide_max <= 1.0 &&
            isfinite(cfg.phase_tracking.cohenite_carbide_max)
        ) || throw(
            ArgumentError(
                "cohenite_carbide_max must be in (0, 1] and finite, got $(cfg.phase_tracking.cohenite_carbide_max)",
            ),
        )
        cfg.phase_tracking.nitride_mode in Set([:roaldite, :carlsbergite, :osbornite]) ||
            throw(
                ArgumentError(
                    "nitride_mode must be :roaldite, :carlsbergite, or :osbornite, got :$(cfg.phase_tracking.nitride_mode)",
                ),
            )
        (
            0.0 < cfg.phase_tracking.r_core_norm < 1.0 &&
            isfinite(cfg.phase_tracking.r_core_norm)
        ) || throw(
            ArgumentError(
                "r_core_norm must be in (0, 1) and finite, got $(cfg.phase_tracking.r_core_norm)",
            ),
        )
        (
            cfg.phase_tracking.r_core_norm < cfg.phase_tracking.r_mantle_norm <= 1.0 &&
            isfinite(cfg.phase_tracking.r_mantle_norm)
        ) || throw(
            ArgumentError(
                "r_mantle_norm must be in (r_core_norm, 1] and finite, got $(cfg.phase_tracking.r_mantle_norm)",
            ),
        )
    end

    if cfg.hydrothermal.active
        (
            0.0 <= cfg.hydrothermal.phi_start < cfg.hydrothermal.phi_end <= 1.0 &&
            isfinite(cfg.hydrothermal.phi_start) &&
            isfinite(cfg.hydrothermal.phi_end)
        ) || throw(
            ArgumentError(
                "phi_start and phi_end must satisfy 0.0 <= phi_start < phi_end <= 1.0 and be finite, got ($(cfg.hydrothermal.phi_start), $(cfg.hydrothermal.phi_end))",
            ),
        )
        @check_positive_finite cfg.hydrothermal.Ra_m_crit
        @check_positive_finite cfg.hydrothermal.Ra_crit
        @check_positive_finite cfg.hydrothermal.c_porous
        @check_positive_finite cfg.hydrothermal.c_free
        @check_positive_finite cfg.hydrothermal.H_layer
        @check_positive_finite cfg.hydrothermal.dT_min
        (
            0.0 < cfg.hydrothermal.k_floor < cfg.hydrothermal.k_cutoff &&
            isfinite(cfg.hydrothermal.k_floor) &&
            isfinite(cfg.hydrothermal.k_cutoff)
        ) || throw(
            ArgumentError(
                "k_floor and k_cutoff must satisfy 0.0 < k_floor < k_cutoff and be finite, got ($(cfg.hydrothermal.k_floor), $(cfg.hydrothermal.k_cutoff))",
            ),
        )
        (
            0.0 < cfg.hydrothermal.picard_damping <= 1.0 &&
            isfinite(cfg.hydrothermal.picard_damping)
        ) || throw(
            ArgumentError(
                "picard_damping must be in (0, 1] and finite, got $(cfg.hydrothermal.picard_damping)",
            ),
        )
        @check_positive_finite cfg.hydrothermal.Pe_crit
        @check_positive_finite cfg.hydrothermal.T_surface_ref
        @check_positive_finite cfg.hydrothermal.gravity
        @check_positive_finite cfg.hydrothermal.cp_fluid
        @check_positive_finite cfg.hydrothermal.alpha_fluid
        @check_positive_finite cfg.hydrothermal.k_fluid_ref
        @check_positive_finite cfg.hydrothermal.rho_fluid_ref
        @check_positive_finite cfg.hydrothermal.mu_fluid_ref
        @check_positive_finite cfg.hydrothermal.kphi_ref
    end

    if cfg.accretion.active
        @check_positive_finite cfg.accretion.M_initial
        @check_positive_finite cfg.accretion.R_initial
        @check_positive_finite cfg.accretion.rho_bulk
        (
            cfg.accretion.M_target >= cfg.accretion.M_initial &&
            isfinite(cfg.accretion.M_target)
        ) || throw(
            ArgumentError(
                "M_target must be >= M_initial and finite, got $(cfg.accretion.M_target)",
            ),
        )
        (
            cfg.accretion.R_target >= cfg.accretion.R_initial &&
            isfinite(cfg.accretion.R_target)
        ) || throw(
            ArgumentError(
                "R_target must be >= R_initial and finite, got $(cfg.accretion.R_target)",
            ),
        )
        @check_nonneg_finite cfg.accretion.t_start_myr
        @check_positive_finite cfg.accretion.t_duration_myr
        @check_positive_finite cfg.accretion.dM_dt_constant
        @check_positive_finite cfg.accretion.dR_dt_constant
        @check_positive_finite cfg.accretion.tau_growth_myr
        r_max_domain = min(
            cfg.geometry.xcenter,
            cfg.geometry.ycenter,
            cfg.grid.xsize - cfg.geometry.xcenter,
            cfg.grid.ysize - cfg.geometry.ycenter,
        )
        r_max_allowed = if cfg.telescoping.active
            (
                isapprox(cfg.geometry.xcenter, cfg.grid.xsize / 2.0; rtol=1e-3) &&
                isapprox(cfg.geometry.ycenter, cfg.grid.ysize / 2.0; rtol=1e-3)
            ) || throw(
                ArgumentError(
                    "Telescoping domain requires planet center (xcenter, ycenter) at domain center (xsize/2, ysize/2)",
                ),
            )
            r_max_domain * (2.0^min(cfg.telescoping.max_telescope_levels, 30))
        else
            r_max_domain
        end
        (cfg.accretion.R_target <= r_max_allowed) || throw(
            ArgumentError(
                "R_target ($(cfg.accretion.R_target)) exceeds distance to domain boundary ($r_max_allowed)",
            ),
        )
        cfg.accretion.mode in Set([
            :constant_rate,
            :linear_radius,
            :exponential,
            :safronov,
            :pebble_bondi,
            :pebble_hill,
            :pebble_auto,
            :multistage,
        ]) || throw(
            ArgumentError(
                "mode must be :constant_rate, :linear_radius, :exponential, :safronov, :pebble_bondi, :pebble_hill, :pebble_auto, or :multistage, got :$(cfg.accretion.mode)",
            ),
        )
        @check_unit_interval cfg.accretion.h_impact
        @check_nonneg_finite cfg.accretion.v_inf
        @check_positive_finite cfg.accretion.cp_rock
        @check_unit_interval cfg.accretion.phi_accreted
        @check_unit_interval cfg.accretion.Xfe_bulk_accreted
        @check_positive_finite cfg.accretion.T_snowline_cond
        @check_unit_interval cfg.accretion.XWsolid_wet
        @check_unit_interval cfg.accretion.XWsolid_dry
        @check_nonneg_finite cfg.accretion.XH2O_wet_wtpct
        @check_nonneg_finite cfg.accretion.XH2O_dry_wtpct
        @check_nonneg_finite cfg.accretion.XC_accreted_ppm
        @check_nonneg_finite cfg.accretion.XN_accreted_ppm
        @check_nonneg_finite cfg.accretion.XS_accreted_ppm
        @check_positive_finite cfg.accretion.stokes_number
        @check_positive_finite cfg.accretion.alpha_turbulence
        @check_nonneg_finite cfg.accretion.Sigma_peb_0
        @check_nonneg_finite cfg.accretion.Sigma_pl_0
        @check_positive_finite cfg.accretion.v_disp_kms
        valid_stage = Set([
            :constant_rate,
            :linear_radius,
            :exponential,
            :safronov,
            :pebble_bondi,
            :pebble_hill,
            :pebble_auto,
        ])
        cfg.accretion.stage1_mode in valid_stage || throw(
            ArgumentError(
                "stage1_mode must be one of $valid_stage, got :$(cfg.accretion.stage1_mode)",
            ),
        )
        cfg.accretion.stage2_mode in valid_stage || throw(
            ArgumentError(
                "stage2_mode must be one of $valid_stage, got :$(cfg.accretion.stage2_mode)",
            ),
        )
        cfg.accretion.stage3_mode in valid_stage || throw(
            ArgumentError(
                "stage3_mode must be one of $valid_stage, got :$(cfg.accretion.stage3_mode)",
            ),
        )
        @check_positive_finite cfg.accretion.f_onset
        @check_positive_finite cfg.accretion.f_iso
        (
            isnan(cfg.accretion.M_onset) ||
            (cfg.accretion.M_onset > 0.0 && isfinite(cfg.accretion.M_onset))
        ) || throw(
            ArgumentError(
                "M_onset must be NaN or > 0 and finite, got $(cfg.accretion.M_onset)"
            ),
        )
        (isnan(cfg.accretion.M_iso) || isfinite(cfg.accretion.M_iso)) ||
            throw(ArgumentError("M_iso must be NaN or finite, got $(cfg.accretion.M_iso)"))
        if !isnan(cfg.accretion.M_onset) &&
            !isnan(cfg.accretion.M_iso) &&
            cfg.accretion.M_iso > 0.0
            cfg.accretion.M_onset < cfg.accretion.M_iso || throw(
                ArgumentError(
                    "M_onset ($(cfg.accretion.M_onset)) must be < M_iso ($(cfg.accretion.M_iso))",
                ),
            )
        end
        (
            0.0 <= cfg.accretion.transition_width <= 0.5 &&
            isfinite(cfg.accretion.transition_width)
        ) || throw(
            ArgumentError(
                "transition_width must be in [0, 0.5] and finite, got $(cfg.accretion.transition_width)",
            ),
        )
    end

    if cfg.telescoping.active
        (
            0.0 < cfg.telescoping.r_threshold_fraction < 1.0 &&
            isfinite(cfg.telescoping.r_threshold_fraction)
        ) || throw(
            ArgumentError(
                "r_threshold_fraction must be in (0, 1) and finite, got $(cfg.telescoping.r_threshold_fraction)",
            ),
        )
        isodd(cfg.grid.Nx) || throw(
            ArgumentError(
                "Telescoping domain requires odd grid.Nx for symmetric centering, got $(cfg.grid.Nx)",
            ),
        )
        isodd(cfg.grid.Ny) || throw(
            ArgumentError(
                "Telescoping domain requires odd grid.Ny for symmetric centering, got $(cfg.grid.Ny)",
            ),
        )
        (
            cfg.telescoping.max_telescope_levels >= 1 &&
            cfg.telescoping.max_telescope_levels <= 30
        ) || throw(
            ArgumentError(
                "max_telescope_levels must be between 1 and 30, got $(cfg.telescoping.max_telescope_levels)",
            ),
        )
        (cfg.telescoping.target_radius > 0.0 && isfinite(cfg.telescoping.target_radius)) ||
            throw(
                ArgumentError(
                    "target_radius must be > 0 and finite, got $(cfg.telescoping.target_radius)",
                ),
            )
        (cfg.telescoping.buffer_markers_per_cell >= 1) || throw(
            ArgumentError(
                "buffer_markers_per_cell must be >= 1, got $(cfg.telescoping.buffer_markers_per_cell)",
            ),
        )
    end

    # Refractory checks
    @check_unit_interval cfg.refractory.f_refr_C
    @check_unit_interval cfg.refractory.f_refr_N
    @check_unit_interval cfg.refractory.f_refr_S
    @check_unit_interval cfg.refractory.f_refr_P
    @check_unit_interval cfg.refractory.f_refr_H
    @check_nonneg_finite cfg.refractory.T_pyrolysis_C
    @check_nonneg_finite cfg.refractory.T_dehydrate_H
    @check_nonneg_finite cfg.refractory.T_pyro_min
    if cfg.refractory.kinetics_active
        @check_positive_finite cfg.refractory.A_C
        @check_positive_finite cfg.refractory.Ea_C
        @check_positive_finite cfg.refractory.A_N
        @check_positive_finite cfg.refractory.Ea_N
        @check_positive_finite cfg.refractory.A_H
        @check_positive_finite cfg.refractory.Ea_H
        @check_nonneg_finite cfg.refractory.dh_pyro_C
        @check_nonneg_finite cfg.refractory.dh_pyro_N
        @check_nonneg_finite cfg.refractory.dh_pyro_H
    end

    # Redox checks
    if cfg.redox.active
        cfg.redox.reference in Set([:mantle, :crust]) || throw(
            ArgumentError(
                "redox.reference must be :mantle or :crust, got $(cfg.redox.reference)"
            ),
        )
        @check_finite cfg.redox.deltaIW_min
        @check_finite cfg.redox.deltaIW_max
        cfg.redox.deltaIW_min <= cfg.redox.deltaIW_max || throw(
            ArgumentError(
                "redox deltaIW_min ($(cfg.redox.deltaIW_min)) must be <= deltaIW_max ($(cfg.redox.deltaIW_max))",
            ),
        )
        @check_unit_interval cfg.redox.initial_x_ferric
    end

    # Volatile mixture checks
    @check_nonneg_finite cfg.volatile_mixture.T_eutectic_ammonia
    @check_nonneg_finite cfg.volatile_mixture.T_freeze_floor
    @check_nonneg_finite cfg.volatile_mixture.alpha_P
    @check_positive_finite cfg.volatile_mixture.P_ref
    @check_unit_interval cfg.volatile_mixture.X_ice_H2O
    @check_unit_interval cfg.volatile_mixture.X_ice_CO2
    @check_unit_interval cfg.volatile_mixture.X_ice_CO
    @check_unit_interval cfg.volatile_mixture.X_ice_CH4
    @check_unit_interval cfg.volatile_mixture.X_ice_NH3
    @check_unit_interval cfg.volatile_mixture.X_ice_N2
    @check_unit_interval cfg.volatile_mixture.X_ice_H2S
    @check_unit_interval cfg.volatile_mixture.X_ice_PH3
    if !(
        cfg.volatile_mixture.T_cond_H2O >= cfg.volatile_mixture.T_cond_NH3 &&
        cfg.volatile_mixture.T_cond_NH3 >= cfg.volatile_mixture.T_cond_CO2 &&
        cfg.volatile_mixture.T_cond_NH3 >= cfg.volatile_mixture.T_cond_H2S &&
        cfg.volatile_mixture.T_cond_CO2 >= cfg.volatile_mixture.T_cond_CH4 &&
        cfg.volatile_mixture.T_cond_H2S >= cfg.volatile_mixture.T_cond_CH4 &&
        cfg.volatile_mixture.T_cond_CH4 >= cfg.volatile_mixture.T_cond_PH3 &&
        cfg.volatile_mixture.T_cond_PH3 >= cfg.volatile_mixture.T_cond_CO &&
        cfg.volatile_mixture.T_cond_CO >= cfg.volatile_mixture.T_cond_N2
    )
        throw(
            ArgumentError(
                "Unphysical snowline ordering: condensation temperatures must satisfy T_cond_H2O >= T_cond_NH3 >= max(T_cond_CO2, T_cond_H2S) >= min(T_cond_CO2, T_cond_H2S) >= T_cond_CH4 >= T_cond_PH3 >= T_cond_CO >= T_cond_N2",
            ),
        )
    end

    # AtmosphereConfig validation
    if cfg.atmosphere.active
        @check_positive_finite cfg.atmosphere.kappa_ir_default
        @check_positive_finite cfg.atmosphere.kappa_vis_default
        cfg.atmosphere.mode in Set([:guillot, :grey, :isothermal]) || throw(
            ArgumentError(
                "atmosphere.mode must be one of :guillot, :grey, :isothermal, got $(cfg.atmosphere.mode)",
            ),
        )
        (0.0 <= cfg.atmosphere.albedo < 1.0 && isfinite(cfg.atmosphere.albedo)) || throw(
            ArgumentError(
                "atmosphere.albedo must be in [0, 1) and finite, got $(cfg.atmosphere.albedo)",
            ),
        )
        @check_positive_finite cfg.atmosphere.gamma_guillot
        @check_positive_finite cfg.atmosphere.T_skin_floor
        (0.0 < cfg.atmosphere.f_rec <= 1.0 && isfinite(cfg.atmosphere.f_rec)) || throw(
            ArgumentError(
                "atmosphere.f_rec must be in (0, 1] and finite, got $(cfg.atmosphere.f_rec)",
            ),
        )
        @check_positive_finite cfg.atmosphere.tau_boil
        @check_positive_finite cfg.atmosphere.b_diff_ref
        for (sp, kap) in cfg.atmosphere.opacities
            (kap >= 0.0 && isfinite(kap)) ||
                throw(ArgumentError("opacity for $sp must be >= 0 and finite, got $kap"))
        end
    end

    # Magma ocean degassing validation
    if cfg.magma_degassing.active
        cfg.magma_degassing.mode in Set([:equilibrium, :dynamic_flux]) || throw(
            ArgumentError(
                "magma_degassing.mode must be :equilibrium or :dynamic_flux, got $(cfg.magma_degassing.mode)",
            ),
        )
        @check_unit_interval cfg.magma_degassing.F_melt_threshold
        @check_unit_interval cfg.magma_degassing.degas_depth_fraction
        (
            0.0 < cfg.magma_degassing.efficiency <= 1.0 &&
            isfinite(cfg.magma_degassing.efficiency)
        ) || throw(
            ArgumentError(
                "magma_degassing.efficiency must be in (0, 1], got $(cfg.magma_degassing.efficiency)",
            ),
        )
    end

    # EscapeConfig validation
    @check_positive_finite cfg.escape.epsilon_xuv
    @check_nonneg_finite cfg.escape.F_xuv_1au_sat
    @check_positive_finite cfg.escape.t_sat_yr
    @check_nonneg_finite cfg.escape.beta_xuv
    (cfg.escape.r_xuv_ratio >= 1.0 && isfinite(cfg.escape.r_xuv_ratio)) || throw(
        ArgumentError(
            "escape.r_xuv_ratio must be >= 1.0 and finite, got $(cfg.escape.r_xuv_ratio)",
        ),
    )

    return nothing
end

"""
Helper function to convert TOML-parsed dictionary into a typed struct with defaults.
"""
function _dict_to_struct(::Type{T}, d::Dict{String,Any}, defaults::T) where {T}
    # Check for unknown / misspelled keys
    for k in keys(d)
        if !hasfield(T, Symbol(k))
            throw(
                ArgumentError("Unknown configuration key '$k' in [$(nameof(T))] section.")
            )
        end
    end

    kwargs = Dict{Symbol,Any}()
    for fname in fieldnames(T)
        sname = String(fname)
        ftype = fieldtype(T, fname)
        if haskey(d, sname)
            val = d[sname]
            if ftype <: SVector
                expected_len = length(ftype)
                if !(val isa AbstractVector) || length(val) != expected_len
                    throw(
                        ArgumentError(
                            "Field '$sname' in [$(nameof(T))] must have exactly $expected_len elements, got $(val)",
                        ),
                    )
                end
                kwargs[fname] = SVector{expected_len,eltype(ftype)}(val)
            elseif ftype <: Real && !(val isa ftype)
                if val isa AbstractString && ftype <: AbstractFloat
                    parsed = tryparse(ftype, val)
                    if parsed !== nothing
                        kwargs[fname] = parsed
                    else
                        throw(
                            ArgumentError(
                                "Cannot parse string '$val' as $ftype for field '$sname'"
                            ),
                        )
                    end
                else
                    kwargs[fname] = convert(ftype, val)
                end
            elseif ftype === Symbol && val isa AbstractString
                kwargs[fname] = Symbol(val)
            elseif ftype === Vector{Symbol} && val isa AbstractVector
                kwargs[fname] = [Symbol(x) for x in val]
            elseif ftype <: AbstractDict && val isa AbstractDict
                K = keytype(ftype)
                V = valtype(ftype)
                kwargs[fname] = Dict{K,V}(
                    (K === Symbol ? Symbol(k) : convert(K, k)) => convert(V, v) for
                    (k, v) in val
                )
            else
                kwargs[fname] = val
            end
        else
            kwargs[fname] = getfield(defaults, fname)
        end
    end
    return T(; kwargs...)
end

const VALID_SECTIONS = Set([
    "grid",
    "geometry",
    "time",
    "solver",
    "poroelasticity",
    "thermodynamics",
    "reaction",
    "materials",
    "output",
    "disk",
    "melting",
    "venting",
    "volatiles",
    "retention",
    "escape",
    "coreformation",
    "metal_partition",
    "phase_tracking",
    "hydrothermal",
    "accretion",
    "telescoping",
    "refractory",
    "volatile_mixture",
    "atmosphere",
    "magma_transport",
    "redox",
    "magma_degassing",
    "mpi",
])

"""
Loads, merges, and validates a `SimulationConfig` from a TOML file or string.

$(SIGNATURES)

# Arguments
- `source`: File path to `.toml` file, or raw TOML string.

# Returns
- `cfg::SimulationConfig`: Validated configuration struct.
"""
function load_config(source::AbstractString)::SimulationConfig
    has_newlines = occursin('\n', source)
    # Resolve file path: direct path or relative to Erebus package root
    resolved_path = if !has_newlines && isfile(source)
        source
    elseif !has_newlines && isfile(joinpath(@__DIR__, "..", source))
        normpath(joinpath(@__DIR__, "..", source))
    else
        nothing
    end

    parsed = if resolved_path !== nothing
        TOML.parsefile(resolved_path)
    elseif endswith(source, ".toml")
        throw(SystemError("opening configuration file: '$source'", 2))
    else
        TOML.parse(source)
    end

    for sec in keys(parsed)
        if sec ∉ VALID_SECTIONS
            throw(ArgumentError("Unknown section '[$sec]' in configuration."))
        end
    end

    def = default_config()

    grid = if haskey(parsed, "grid")
        _dict_to_struct(GridConfig, parsed["grid"], def.grid)
    else
        def.grid
    end
    geom = if haskey(parsed, "geometry")
        _dict_to_struct(GeometryConfig, parsed["geometry"], def.geometry)
    else
        def.geometry
    end
    time = if haskey(parsed, "time")
        _dict_to_struct(TimeConfig, parsed["time"], def.time)
    else
        def.time
    end
    solv = if haskey(parsed, "solver")
        _dict_to_struct(SolverConfig, parsed["solver"], def.solver)
    else
        def.solver
    end
    poro = if haskey(parsed, "poroelasticity")
        _dict_to_struct(PoroelasticConfig, parsed["poroelasticity"], def.poroelasticity)
    else
        def.poroelasticity
    end
    therm = if haskey(parsed, "thermodynamics")
        _dict_to_struct(ThermalConfig, parsed["thermodynamics"], def.thermodynamics)
    else
        def.thermodynamics
    end
    react = if haskey(parsed, "reaction")
        _dict_to_struct(ReactionConfig, parsed["reaction"], def.reaction)
    else
        def.reaction
    end
    mat = if haskey(parsed, "materials")
        _dict_to_struct(MaterialConfig, parsed["materials"], def.materials)
    else
        def.materials
    end
    out = if haskey(parsed, "output")
        _dict_to_struct(OutputConfig, parsed["output"], def.output)
    else
        def.output
    end
    dsk = if haskey(parsed, "disk")
        _dict_to_struct(DiskConfig, parsed["disk"], def.disk)
    else
        def.disk
    end
    melt = if haskey(parsed, "melting")
        _dict_to_struct(MeltingConfig, parsed["melting"], def.melting)
    else
        def.melting
    end
    vent = if haskey(parsed, "venting")
        _dict_to_struct(VentingConfig, parsed["venting"], def.venting)
    else
        def.venting
    end
    vol = if haskey(parsed, "volatiles")
        _dict_to_struct(VolatilesConfig, parsed["volatiles"], def.volatiles)
    else
        def.volatiles
    end
    ret = if haskey(parsed, "retention")
        _dict_to_struct(RetentionConfig, parsed["retention"], def.retention)
    else
        def.retention
    end
    esc = if haskey(parsed, "escape")
        parsed_esc = parsed["escape"]
        def_esc =
            if !haskey(parsed_esc, "R_planet") && geom.rplanet != def.geometry.rplanet
                EscapeConfig(;
                    active=def.escape.active,
                    M_planet=def.escape.M_planet,
                    R_planet=geom.rplanet,
                    T_exobase=def.escape.T_exobase,
                    R_exobase=geom.rplanet,
                    species=def.escape.species,
                    multi_species=def.escape.multi_species,
                    species_list=def.escape.species_list,
                    gamma=def.escape.gamma,
                    hydrodynamic=def.escape.hydrodynamic,
                    xuv_driven=def.escape.xuv_driven,
                    epsilon_xuv=def.escape.epsilon_xuv,
                    F_xuv_1au_sat=def.escape.F_xuv_1au_sat,
                    t_sat_yr=def.escape.t_sat_yr,
                    beta_xuv=def.escape.beta_xuv,
                    r_xuv_ratio=def.escape.r_xuv_ratio,
                    tidal_correction=def.escape.tidal_correction,
                )
            else
                def.escape
            end
        _dict_to_struct(EscapeConfig, parsed_esc, def_esc)
    else
        if geom.rplanet != def.geometry.rplanet
            EscapeConfig(;
                active=def.escape.active,
                M_planet=def.escape.M_planet,
                R_planet=geom.rplanet,
                T_exobase=def.escape.T_exobase,
                R_exobase=geom.rplanet,
                species=def.escape.species,
                multi_species=def.escape.multi_species,
                species_list=def.escape.species_list,
                gamma=def.escape.gamma,
                hydrodynamic=def.escape.hydrodynamic,
                xuv_driven=def.escape.xuv_driven,
                epsilon_xuv=def.escape.epsilon_xuv,
                F_xuv_1au_sat=def.escape.F_xuv_1au_sat,
                t_sat_yr=def.escape.t_sat_yr,
                beta_xuv=def.escape.beta_xuv,
                r_xuv_ratio=def.escape.r_xuv_ratio,
                tidal_correction=def.escape.tidal_correction,
            )
        else
            def.escape
        end
    end
    coreform = if haskey(parsed, "coreformation")
        _dict_to_struct(CoreFormationConfig, parsed["coreformation"], def.coreformation)
    else
        def.coreformation
    end
    metal_part = if haskey(parsed, "metal_partition")
        _dict_to_struct(MetalPartitionConfig, parsed["metal_partition"], def.metal_partition)
    else
        def.metal_partition
    end
    phase_track = if haskey(parsed, "phase_tracking")
        _dict_to_struct(PhaseTrackingConfig, parsed["phase_tracking"], def.phase_tracking)
    else
        def.phase_tracking
    end
    hydrotherm = if haskey(parsed, "hydrothermal")
        _dict_to_struct(HydrothermalConfig, parsed["hydrothermal"], def.hydrothermal)
    else
        def.hydrothermal
    end
    acc = if haskey(parsed, "accretion")
        _dict_to_struct(AccretionConfig, parsed["accretion"], def.accretion)
    else
        def.accretion
    end
    tele = if haskey(parsed, "telescoping")
        _dict_to_struct(TelescopingConfig, parsed["telescoping"], def.telescoping)
    else
        def.telescoping
    end
    refr = if haskey(parsed, "refractory")
        _dict_to_struct(RefractoryConfig, parsed["refractory"], def.refractory)
    else
        def.refractory
    end
    volmix = if haskey(parsed, "volatile_mixture")
        _dict_to_struct(VolatileMixtureConfig, parsed["volatile_mixture"], def.volatile_mixture)
    else
        def.volatile_mixture
    end
    atm = if haskey(parsed, "atmosphere")
        _dict_to_struct(AtmosphereConfig, parsed["atmosphere"], def.atmosphere)
    else
        def.atmosphere
    end
    magma = if haskey(parsed, "magma_transport")
        _dict_to_struct(MagmaTransportConfig, parsed["magma_transport"], def.magma_transport)
    else
        def.magma_transport
    end
    rdx = if haskey(parsed, "redox")
        _dict_to_struct(RedoxConfig, parsed["redox"], def.redox)
    else
        def.redox
    end
    magma_degas = if haskey(parsed, "magma_degassing")
        _dict_to_struct(
            MagmaOceanDegassingConfig, parsed["magma_degassing"], def.magma_degassing
        )
    else
        def.magma_degassing
    end
    mpi = if haskey(parsed, "mpi")
        _dict_to_struct(MPIConfig, parsed["mpi"], def.mpi)
    else
        def.mpi
    end

    cfg = SimulationConfig(;
        grid=grid,
        geometry=geom,
        time=time,
        solver=solv,
        poroelasticity=poro,
        thermodynamics=therm,
        reaction=react,
        materials=mat,
        output=out,
        disk=dsk,
        melting=melt,
        venting=vent,
        volatiles=vol,
        retention=ret,
        escape=esc,
        coreformation=coreform,
        metal_partition=metal_part,
        phase_tracking=phase_track,
        hydrothermal=hydrotherm,
        accretion=acc,
        telescoping=tele,
        refractory=refr,
        volatile_mixture=volmix,
        atmosphere=atm,
        magma_transport=magma,
        redox=rdx,
        magma_degassing=magma_degas,
        mpi=mpi,
    )

    validate_config(cfg)
    return cfg
end

"""
Converts a struct into a Dict suitable for TOML serialization.
"""
function _struct_to_dict(s)
    d = Dict{String,Any}()
    for fname in fieldnames(typeof(s))
        val = getfield(s, fname)
        if val isa SVector
            d[String(fname)] = collect(val)
        elseif val isa Symbol
            d[String(fname)] = String(val)
        elseif val isa AbstractVector{Symbol}
            d[String(fname)] = [String(x) for x in val]
        elseif val isa AbstractDict{Symbol,<:Any}
            d[String(fname)] = Dict{String,Any}(String(k) => v for (k, v) in val)
        else
            d[String(fname)] = val
        end
    end
    return d
end

"""
Converts a `SimulationConfig` to a nested dictionary representation.

$(SIGNATURES)
"""
function config_to_dict(cfg::SimulationConfig)::Dict{String,Any}
    return Dict{String,Any}(
        "grid" => _struct_to_dict(cfg.grid),
        "geometry" => _struct_to_dict(cfg.geometry),
        "time" => _struct_to_dict(cfg.time),
        "solver" => _struct_to_dict(cfg.solver),
        "poroelasticity" => _struct_to_dict(cfg.poroelasticity),
        "thermodynamics" => _struct_to_dict(cfg.thermodynamics),
        "reaction" => _struct_to_dict(cfg.reaction),
        "materials" => _struct_to_dict(cfg.materials),
        "output" => _struct_to_dict(cfg.output),
        "disk" => _struct_to_dict(cfg.disk),
        "melting" => _struct_to_dict(cfg.melting),
        "venting" => _struct_to_dict(cfg.venting),
        "volatiles" => _struct_to_dict(cfg.volatiles),
        "retention" => _struct_to_dict(cfg.retention),
        "escape" => _struct_to_dict(cfg.escape),
        "coreformation" => _struct_to_dict(cfg.coreformation),
        "metal_partition" => _struct_to_dict(cfg.metal_partition),
        "phase_tracking" => _struct_to_dict(cfg.phase_tracking),
        "hydrothermal" => _struct_to_dict(cfg.hydrothermal),
        "accretion" => _struct_to_dict(cfg.accretion),
        "telescoping" => _struct_to_dict(cfg.telescoping),
        "refractory" => _struct_to_dict(cfg.refractory),
        "volatile_mixture" => _struct_to_dict(cfg.volatile_mixture),
        "atmosphere" => _struct_to_dict(cfg.atmosphere),
        "magma_transport" => _struct_to_dict(cfg.magma_transport),
        "redox" => _struct_to_dict(cfg.redox),
        "magma_degassing" => _struct_to_dict(cfg.magma_degassing),
        "mpi" => _struct_to_dict(cfg.mpi),
    )
end

"""
Saves a `SimulationConfig` to an IO stream or a `.toml` file.

$(SIGNATURES)

# Arguments
- `io_or_path`: Output IO stream or file path.
- `cfg`: Configuration to serialize.
"""
function save_config(io::IO, cfg::SimulationConfig)
    d = config_to_dict(cfg)
    TOML.print(io, d; sorted=true)
    return io
end

function save_config(path::AbstractString, cfg::SimulationConfig)
    open(path, "w") do io
        return save_config(io, cfg)
    end
    return path
end

function save_config(cfg::SimulationConfig)::String
    io = IOBuffer()
    save_config(io, cfg)
    return String(take!(io))
end

"""
Serializes a `SimulationConfig` to a TOML-formatted string.

$(SIGNATURES)
"""
function serialize_config(cfg::SimulationConfig)::String
    return save_config(cfg)
end

"""
Parses a `SimulationConfig` from a TOML-formatted string.

$(SIGNATURES)
"""
function parse_config_string(str::AbstractString)::SimulationConfig
    return load_config(str)
end
