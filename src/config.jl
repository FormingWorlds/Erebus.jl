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
Planetesimal surface volatile degassing and venting parameters.

$(FIELDS)
"""
Base.@kwdef struct VentingConfig
    active::Bool = false
    mode::Symbol = :darcy_sink
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

Configures the standalone thermodynamic speciation and volatile solubility library
(`Erebus.Physics.Volatiles`). Dynamic reactive transport coupling in 2D fluid flow is in development.

$(FIELDS)
"""
Base.@kwdef struct VolatilesConfig
    active::Bool = false
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
Atmospheric Jeans kinetic escape and volatile mass loss parameters.

$(FIELDS)
"""
Base.@kwdef struct EscapeConfig
    active::Bool = false
    M_planet::Float64 = 1.309e18
    R_planet::Float64 = 50_000.0
    T_exobase::Float64 = 200.0
    R_exobase::Float64 = 50_000.0
    species::Symbol = :H2O
end

"""
Iron core formation by metal percolation and gravitational settling parameters.

$(FIELDS)
"""
Base.@kwdef struct CoreFormationConfig
    percolation_active::Bool = false
    settling_active::Bool = false
    rho_metal::Float64 = 7200.0
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
    droplet_size_mode::Symbol = :weber_mean
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
    escape::EscapeConfig = EscapeConfig()
    coreformation::CoreFormationConfig = CoreFormationConfig()
end

"""
Returns the default simulation configuration matching baseline constants.

$(SIGNATURES)
"""
function default_config()::SimulationConfig
    return SimulationConfig()
end

"""
Validates physical bounds and numerical consistency of a `SimulationConfig`.

$(SIGNATURES)

# Raises
- `ArgumentError` if any configuration parameter violates physical bounds or compiled grid constraints.
"""
function validate_config(cfg::SimulationConfig)
    # Grid checks: must be valid dimensions and domain bounds
    cfg.grid.Nx >= 3 || throw(ArgumentError("Grid Nx must be >= 3, got $(cfg.grid.Nx)"))
    cfg.grid.Ny >= 3 || throw(ArgumentError("Grid Ny must be >= 3, got $(cfg.grid.Ny)"))
    cfg.grid.xsize > 0.0 ||
        throw(ArgumentError("Domain xsize must be > 0, got $(cfg.grid.xsize)"))
    cfg.grid.ysize > 0.0 ||
        throw(ArgumentError("Domain ysize must be > 0, got $(cfg.grid.ysize)"))

    # Geometry checks
    cfg.geometry.rplanet > 0.0 ||
        throw(ArgumentError("Planet radius must be > 0, got $(cfg.geometry.rplanet)"))
    cfg.geometry.rcrust > 0.0 ||
        throw(ArgumentError("Crust radius must be > 0, got $(cfg.geometry.rcrust)"))
    cfg.geometry.rcrust <= cfg.geometry.rplanet ||
        throw(ArgumentError("Crust radius must be <= planet radius"))
    cfg.geometry.metric_regularization_cells > 0.0 &&
    isfinite(cfg.geometry.metric_regularization_cells) || throw(
        ArgumentError(
            "metric_regularization_cells must be > 0 and finite, got $(cfg.geometry.metric_regularization_cells)",
        ),
    )
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
    cfg.time.dt_initial > 0.0 ||
        throw(ArgumentError("Initial dt must be > 0, got $(cfg.time.dt_initial)"))
    cfg.time.dt_longest >= cfg.time.dt_initial ||
        throw(ArgumentError("dt_longest must be >= dt_initial"))
    cfg.time.n_steps >= 1 ||
        throw(ArgumentError("n_steps must be >= 1, got $(cfg.time.n_steps)"))
    cfg.time.start_time >= 0.0 ||
        throw(ArgumentError("start_time must be >= 0, got $(cfg.time.start_time)"))
    cfg.time.endtime > cfg.time.start_time ||
        throw(ArgumentError("endtime must be > start_time"))
    cfg.time.start_step >= 1 ||
        throw(ArgumentError("start_step must be >= 1, got $(cfg.time.start_step)"))
    isfinite(cfg.time.dt_initial) || throw(ArgumentError("dt_initial must be finite"))
    isfinite(cfg.time.dt_longest) || throw(ArgumentError("dt_longest must be finite"))

    # Poroelasticity checks
    cfg.poroelasticity.betasolid >= 0.0 ||
        throw(ArgumentError("betasolid must be >= 0, got $(cfg.poroelasticity.betasolid)"))
    cfg.poroelasticity.betafluid >= 0.0 ||
        throw(ArgumentError("betafluid must be >= 0, got $(cfg.poroelasticity.betafluid)"))
    isfinite(cfg.poroelasticity.betasolid) ||
        throw(ArgumentError("betasolid must be finite"))
    isfinite(cfg.poroelasticity.betafluid) ||
        throw(ArgumentError("betafluid must be finite"))
    0.0 < cfg.poroelasticity.phimin < cfg.poroelasticity.phimax < 1.0 || throw(
        ArgumentError(
            "Porosity bounds must satisfy 0 < phimin < phimax < 1, got phimin=$(cfg.poroelasticity.phimin), phimax=$(cfg.poroelasticity.phimax)",
        ),
    )
    cfg.poroelasticity.kappa_frac >= 0.0 && isfinite(cfg.poroelasticity.kappa_frac) ||
        throw(
            ArgumentError(
                "kappa_frac must be >= 0 and finite, got $(cfg.poroelasticity.kappa_frac)"
            ),
        )
    cfg.poroelasticity.gamma_frac > 0.0 && isfinite(cfg.poroelasticity.gamma_frac) || throw(
        ArgumentError(
            "gamma_frac must be > 0 and finite, got $(cfg.poroelasticity.gamma_frac)"
        ),
    )
    cfg.poroelasticity.k_frac_max > 0.0 && isfinite(cfg.poroelasticity.k_frac_max) || throw(
        ArgumentError(
            "k_frac_max must be > 0 and finite, got $(cfg.poroelasticity.k_frac_max)"
        ),
    )

    # Solver checks
    cfg.solver.titermax >= 1 ||
        throw(ArgumentError("titermax must be >= 1, got $(cfg.solver.titermax)"))
    cfg.solver.nplast >= 1 ||
        throw(ArgumentError("nplast must be >= 1, got $(cfg.solver.nplast)"))
    cfg.solver.titermax <= cfg.solver.nplast || throw(
        ArgumentError(
            "titermax ($(cfg.solver.titermax)) must be <= nplast ($(cfg.solver.nplast)) to prevent array bounds overflow in plastic convergence tracking",
        ),
    )
    cfg.solver.etamin > 0.0 ||
        throw(ArgumentError("etamin must be > 0, got $(cfg.solver.etamin)"))
    cfg.solver.etamax >= cfg.solver.etamin ||
        throw(ArgumentError("etamax must be >= etamin"))
    cfg.solver.etaphikoef > 0.0 ||
        throw(ArgumentError("etaphikoef must be > 0, got $(cfg.solver.etaphikoef)"))

    # Output checks
    cfg.output.savematstep >= 1 ||
        throw(ArgumentError("savematstep must be >= 1, got $(cfg.output.savematstep)"))
    cfg.output.visstep >= 1 ||
        throw(ArgumentError("visstep must be >= 1, got $(cfg.output.visstep)"))
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
    0.0 <= cfg.thermodynamics.ratio_al <= 1.0 || throw(
        ArgumentError("ratio_al must be in [0, 1], got $(cfg.thermodynamics.ratio_al)")
    )
    0.0 <= cfg.thermodynamics.ratio_fe <= 1.0 || throw(
        ArgumentError("ratio_fe must be in [0, 1], got $(cfg.thermodynamics.ratio_fe)")
    )
    cfg.thermodynamics.tmfluidphase < cfg.thermodynamics.tmsolidphase || throw(
        ArgumentError(
            "tmfluidphase ($(cfg.thermodynamics.tmfluidphase)) must be < tmsolidphase ($(cfg.thermodynamics.tmsolidphase))",
        ),
    )
    cfg.thermodynamics.Lᶠ > 0.0 ||
        throw(ArgumentError("Lᶠ must be > 0, got $(cfg.thermodynamics.Lᶠ)"))
    cfg.thermodynamics.E_al > 0.0 && isfinite(cfg.thermodynamics.E_al) ||
        throw(ArgumentError("E_al must be > 0 and finite"))
    cfg.thermodynamics.f_al > 0.0 && isfinite(cfg.thermodynamics.f_al) ||
        throw(ArgumentError("f_al must be > 0 and finite"))
    cfg.thermodynamics.t_half_al > 0.0 && isfinite(cfg.thermodynamics.t_half_al) ||
        throw(ArgumentError("t_half_al must be > 0 and finite"))
    cfg.thermodynamics.E_fe > 0.0 && isfinite(cfg.thermodynamics.E_fe) ||
        throw(ArgumentError("E_fe must be > 0 and finite"))
    cfg.thermodynamics.f_fe > 0.0 && isfinite(cfg.thermodynamics.f_fe) ||
        throw(ArgumentError("f_fe must be > 0 and finite"))
    cfg.thermodynamics.t_half_fe > 0.0 && isfinite(cfg.thermodynamics.t_half_fe) ||
        throw(ArgumentError("t_half_fe must be > 0 and finite"))
    cfg.thermodynamics.fluid_viscosity_mode in Set([:arrhenius, :constant]) || throw(
        ArgumentError(
            "fluid_viscosity_mode must be :arrhenius or :constant, got $(cfg.thermodynamics.fluid_viscosity_mode)",
        ),
    )
    cfg.thermodynamics.fluid_viscosity_Ea >= 0.0 &&
    isfinite(cfg.thermodynamics.fluid_viscosity_Ea) ||
        throw(ArgumentError("fluid_viscosity_Ea must be >= 0 and finite"))
    cfg.thermodynamics.fluid_viscosity_T0 > 0.0 &&
    isfinite(cfg.thermodynamics.fluid_viscosity_T0) ||
        throw(ArgumentError("fluid_viscosity_T0 must be > 0 and finite"))
    cfg.thermodynamics.fluid_viscosity_eta0 > 0.0 &&
    isfinite(cfg.thermodynamics.fluid_viscosity_eta0) ||
        throw(ArgumentError("fluid_viscosity_eta0 must be > 0 and finite"))
    0.0 <= cfg.thermodynamics.emissivity <= 1.0 || throw(
        ArgumentError("emissivity must be in [0, 1], got $(cfg.thermodynamics.emissivity)"),
    )
    cfg.thermodynamics.sigma_sb > 0.0 && isfinite(cfg.thermodynamics.sigma_sb) || throw(
        ArgumentError(
            "sigma_sb must be > 0 and finite, got $(cfg.thermodynamics.sigma_sb)"
        ),
    )

    # Disk checks
    cfg.disk.model in Set([:fixed, :monotonic, :class1_to_class2, :class0_to_class2]) ||
        throw(
            ArgumentError(
                "disk model must be :fixed, :monotonic, :class1_to_class2, or :class0_to_class2, got $(cfg.disk.model)",
            ),
        )
    cfg.disk.t_ambient > 0.0 && isfinite(cfg.disk.t_ambient) ||
        throw(ArgumentError("t_ambient must be > 0 and finite, got $(cfg.disk.t_ambient)"))
    cfg.disk.orbital_distance_au > 0.0 && isfinite(cfg.disk.orbital_distance_au) || throw(
        ArgumentError(
            "orbital_distance_au must be > 0 and finite, got $(cfg.disk.orbital_distance_au)",
        ),
    )
    cfg.disk.stellar_mass_msun > 0.0 && isfinite(cfg.disk.stellar_mass_msun) || throw(
        ArgumentError(
            "stellar_mass_msun must be > 0 and finite, got $(cfg.disk.stellar_mass_msun)",
        ),
    )
    cfg.disk.t_cloud > 0.0 && isfinite(cfg.disk.t_cloud) ||
        throw(ArgumentError("t_cloud must be > 0 and finite, got $(cfg.disk.t_cloud)"))
    cfg.disk.t_irr_1au > 0.0 && isfinite(cfg.disk.t_irr_1au) ||
        throw(ArgumentError("t_irr_1au must be > 0 and finite, got $(cfg.disk.t_irr_1au)"))
    cfg.disk.t_peak_1au > 0.0 && isfinite(cfg.disk.t_peak_1au) || throw(
        ArgumentError("t_peak_1au must be > 0 and finite, got $(cfg.disk.t_peak_1au)")
    )
    cfg.disk.t_peak_time_1au_myr > 0.0 && isfinite(cfg.disk.t_peak_time_1au_myr) || throw(
        ArgumentError(
            "t_peak_time_1au_myr must be > 0 and finite, got $(cfg.disk.t_peak_time_1au_myr)",
        ),
    )
    cfg.disk.t_visc_0_myr > 0.0 && isfinite(cfg.disk.t_visc_0_myr) || throw(
        ArgumentError("t_visc_0_myr must be > 0 and finite, got $(cfg.disk.t_visc_0_myr)"),
    )
    cfg.disk.gamma > 0.0 && isfinite(cfg.disk.gamma) ||
        throw(ArgumentError("gamma must be > 0 and finite, got $(cfg.disk.gamma)"))
    cfg.disk.alpha > 0.0 && isfinite(cfg.disk.alpha) ||
        throw(ArgumentError("alpha must be > 0 and finite, got $(cfg.disk.alpha)"))
    cfg.disk.q_irr > 0.0 && isfinite(cfg.disk.q_irr) ||
        throw(ArgumentError("q_irr must be > 0 and finite, got $(cfg.disk.q_irr)"))
    cfg.disk.q_visc > 0.0 && isfinite(cfg.disk.q_visc) ||
        throw(ArgumentError("q_visc must be > 0 and finite, got $(cfg.disk.q_visc)"))
    cfg.disk.p_r_t >= 0.0 && isfinite(cfg.disk.p_r_t) ||
        throw(ArgumentError("p_r_t must be >= 0 and finite, got $(cfg.disk.p_r_t)"))
    cfg.disk.p_m_irr >= 0.0 && isfinite(cfg.disk.p_m_irr) ||
        throw(ArgumentError("p_m_irr must be >= 0 and finite, got $(cfg.disk.p_m_irr)"))
    cfg.disk.p_m_visc >= 0.0 && isfinite(cfg.disk.p_m_visc) ||
        throw(ArgumentError("p_m_visc must be >= 0 and finite, got $(cfg.disk.p_m_visc)"))
    cfg.disk.p_m_t >= 0.0 && isfinite(cfg.disk.p_m_t) ||
        throw(ArgumentError("p_m_t must be >= 0 and finite, got $(cfg.disk.p_m_t)"))
    cfg.disk.p_m_visc_decay >= 0.0 && isfinite(cfg.disk.p_m_visc_decay) || throw(
        ArgumentError(
            "p_m_visc_decay must be >= 0 and finite, got $(cfg.disk.p_m_visc_decay)"
        ),
    )
    cfg.disk.t_dispersal_myr > 0.0 && isfinite(cfg.disk.t_dispersal_myr) || throw(
        ArgumentError(
            "t_dispersal_myr must be > 0 and finite, got $(cfg.disk.t_dispersal_myr)"
        ),
    )
    cfg.disk.dt_dispersal_myr > 0.0 && isfinite(cfg.disk.dt_dispersal_myr) || throw(
        ArgumentError(
            "dt_dispersal_myr must be > 0 and finite, got $(cfg.disk.dt_dispersal_myr)"
        ),
    )
    cfg.disk.p_amb_disk > 0.0 && isfinite(cfg.disk.p_amb_disk) || throw(
        ArgumentError("p_amb_disk must be > 0 and finite, got $(cfg.disk.p_amb_disk)")
    )
    cfg.disk.p_amb_space >= 0.0 && isfinite(cfg.disk.p_amb_space) || throw(
        ArgumentError("p_amb_space must be >= 0 and finite, got $(cfg.disk.p_amb_space)"),
    )
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
    cfg.reaction.dtreaction_hydration > 0.0 &&
    isfinite(cfg.reaction.dtreaction_hydration) ||
        throw(ArgumentError("dtreaction_hydration must be > 0 and finite"))
    cfg.reaction.dtreaction_dehydration > 0.0 &&
    isfinite(cfg.reaction.dtreaction_dehydration) ||
        throw(ArgumentError("dtreaction_dehydration must be > 0 and finite"))
    cfg.reaction.delta_H > 0.0 && isfinite(cfg.reaction.delta_H) ||
        throw(ArgumentError("delta_H must be > 0 and finite"))
    cfg.reaction.delta_S > 0.0 && isfinite(cfg.reaction.delta_S) ||
        throw(ArgumentError("delta_S must be > 0 and finite"))
    cfg.reaction.A_I > 0.0 && isfinite(cfg.reaction.A_I) ||
        throw(ArgumentError("A_I must be > 0 and finite"))
    cfg.reaction.b_I > 0.0 && isfinite(cfg.reaction.b_I) ||
        throw(ArgumentError("b_I must be > 0 and finite"))
    cfg.reaction.c_I > 0.0 && isfinite(cfg.reaction.c_I) ||
        throw(ArgumentError("c_I must be > 0 and finite"))
    cfg.reaction.Sxo_B > 0.0 && isfinite(cfg.reaction.Sxo_B) ||
        throw(ArgumentError("Sxo_B must be > 0 and finite"))
    cfg.reaction.Tscl_B > 0.0 && isfinite(cfg.reaction.Tscl_B) ||
        throw(ArgumentError("Tscl_B must be > 0 and finite"))
    cfg.reaction.To_B > 0.0 && isfinite(cfg.reaction.To_B) ||
        throw(ArgumentError("To_B must be > 0 and finite"))
    0.0 < cfg.reaction.alpha_relaxation <= 1.0 || throw(
        ArgumentError(
            "alpha_relaxation must be in (0, 1], got $(cfg.reaction.alpha_relaxation)"
        ),
    )
    0.0 <= cfg.reaction.pfcoeff <= 1.0 ||
        throw(ArgumentError("pfcoeff must be in [0, 1], got $(cfg.reaction.pfcoeff)"))
    cfg.reaction.pferrmax > 0.0 && isfinite(cfg.reaction.pferrmax) ||
        throw(ArgumentError("pferrmax must be > 0 and finite"))
    cfg.reaction.p_cavitation >= 0.0 && isfinite(cfg.reaction.p_cavitation) ||
        throw(ArgumentError("p_cavitation must be >= 0 and finite"))

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
        (cfg.melting.L_melt > 0.0 && isfinite(cfg.melting.L_melt)) || throw(
            ArgumentError(
                "Melting L_melt must be > 0 and finite, got $(cfg.melting.L_melt)"
            ),
        )
        (cfg.melting.rho_melt > 0.0 && isfinite(cfg.melting.rho_melt)) || throw(
            ArgumentError(
                "Melting rho_melt must be > 0 and finite, got $(cfg.melting.rho_melt)"
            ),
        )
        (cfg.melting.alpha_eta >= 0.0 && isfinite(cfg.melting.alpha_eta)) || throw(
            ArgumentError(
                "Melting alpha_eta must be >= 0 and finite, got $(cfg.melting.alpha_eta)",
            ),
        )
        (0.0 < cfg.melting.phi_crit < 1.0 && isfinite(cfg.melting.phi_crit)) || throw(
            ArgumentError(
                "Melting phi_crit must be in (0, 1), got $(cfg.melting.phi_crit)"
            ),
        )
        (cfg.melting.eta_melt > 0.0 && isfinite(cfg.melting.eta_melt)) || throw(
            ArgumentError(
                "Melting eta_melt must be > 0 and finite, got $(cfg.melting.eta_melt)"
            ),
        )
        (cfg.melting.dpdt_clapeyron >= 0.0 && isfinite(cfg.melting.dpdt_clapeyron)) ||
            throw(
                ArgumentError(
                    "Melting dpdt_clapeyron must be >= 0 and finite, got $(cfg.melting.dpdt_clapeyron)",
                ),
            )
        cfg.melting.latent_heat_mode == :apparent_cp || throw(
            ArgumentError(
                "Melting latent_heat_mode must be :apparent_cp, got $(cfg.melting.latent_heat_mode)",
            ),
        )

        if cfg.melting.soft_turbulence
            (cfg.melting.turb_exponent > 0.0 && isfinite(cfg.melting.turb_exponent)) ||
                throw(
                    ArgumentError(
                        "Melting turb_exponent must be > 0 and finite, got $(cfg.melting.turb_exponent)",
                    ),
                )
            (
                cfg.melting.eta_fluid_silicate > 0.0 &&
                isfinite(cfg.melting.eta_fluid_silicate)
            ) || throw(
                ArgumentError(
                    "Melting eta_fluid_silicate must be > 0 and finite, got $(cfg.melting.eta_fluid_silicate)",
                ),
            )
            (0.0 <= cfg.melting.F_turb_start < cfg.melting.F_turb_end <= 1.0) || throw(
                ArgumentError(
                    "Melting F_turb bounds must satisfy 0 <= F_turb_start < F_turb_end <= 1, got [$(cfg.melting.F_turb_start), $(cfg.melting.F_turb_end)]",
                ),
            )
            (cfg.melting.dT_turb_min > 0.0 && isfinite(cfg.melting.dT_turb_min)) || throw(
                ArgumentError(
                    "Melting dT_turb_min must be > 0 and finite, got $(cfg.melting.dT_turb_min)",
                ),
            )
            (cfg.melting.T_surface_ref > 0.0 && isfinite(cfg.melting.T_surface_ref)) ||
                throw(
                    ArgumentError(
                        "Melting T_surface_ref must be > 0 and finite, got $(cfg.melting.T_surface_ref)",
                    ),
                )
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

    # Venting checks
    cfg.venting.mode in Set([:darcy_sink, :hydrofracture_gated]) || throw(
        ArgumentError(
            "venting mode must be :darcy_sink or :hydrofracture_gated, got $(cfg.venting.mode)",
        ),
    )
    cfg.venting.k_vent > 0.0 && isfinite(cfg.venting.k_vent) ||
        throw(ArgumentError("k_vent must be > 0 and finite, got $(cfg.venting.k_vent)"))
    cfg.venting.conductance_factor > 0.0 && isfinite(cfg.venting.conductance_factor) ||
        throw(
            ArgumentError(
                "conductance_factor must be > 0 and finite, got $(cfg.venting.conductance_factor)",
            ),
        )
    cfg.venting.L_sublimation > 0.0 && isfinite(cfg.venting.L_sublimation) || throw(
        ArgumentError(
            "L_sublimation must be > 0 and finite, got $(cfg.venting.L_sublimation)"
        ),
    )
    cfg.venting.t_freeze > 0.0 && isfinite(cfg.venting.t_freeze) ||
        throw(ArgumentError("t_freeze must be > 0 and finite, got $(cfg.venting.t_freeze)"))
    cfg.venting.dt_seal > 0.0 && isfinite(cfg.venting.dt_seal) ||
        throw(ArgumentError("dt_seal must be > 0 and finite, got $(cfg.venting.dt_seal)"))
    (0.0 < cfg.venting.k_seal_min_ratio <= 1.0 && isfinite(cfg.venting.k_seal_min_ratio)) ||
        throw(
            ArgumentError(
                "k_seal_min_ratio must be in (0, 1], got $(cfg.venting.k_seal_min_ratio)"
            ),
        )

    # Volatiles checks
    (isfinite(cfg.volatiles.fO2_delta_IW) && abs(cfg.volatiles.fO2_delta_IW) <= 50.0) ||
        throw(
            ArgumentError(
                "fO2_delta_IW must be finite and within [-50, 50], got $(cfg.volatiles.fO2_delta_IW)",
            ),
        )
    (
        cfg.volatiles.water_solubility_coeff > 0.0 &&
        isfinite(cfg.volatiles.water_solubility_coeff)
    ) || throw(
        ArgumentError(
            "water_solubility_coeff must be > 0 and finite, got $(cfg.volatiles.water_solubility_coeff)",
        ),
    )
    (
        cfg.volatiles.nitrogen_henry_coeff > 0.0 &&
        isfinite(cfg.volatiles.nitrogen_henry_coeff)
    ) || throw(
        ArgumentError(
            "nitrogen_henry_coeff must be > 0 and finite, got $(cfg.volatiles.nitrogen_henry_coeff)",
        ),
    )
    (
        cfg.volatiles.nitrogen_nitride_capacity > 0.0 &&
        isfinite(cfg.volatiles.nitrogen_nitride_capacity)
    ) || throw(
        ArgumentError(
            "nitrogen_nitride_capacity must be > 0 and finite, got $(cfg.volatiles.nitrogen_nitride_capacity)",
        ),
    )
    (cfg.volatiles.t_organic_devol > 0.0 && isfinite(cfg.volatiles.t_organic_devol)) ||
        throw(
            ArgumentError(
                "t_organic_devol must be > 0 and finite, got $(cfg.volatiles.t_organic_devol)",
            ),
        )
    (cfg.volatiles.dt_organic_devol > 0.0 && isfinite(cfg.volatiles.dt_organic_devol)) ||
        throw(
            ArgumentError(
                "dt_organic_devol must be > 0 and finite, got $(cfg.volatiles.dt_organic_devol)",
            ),
        )
    (
        cfg.volatiles.organic_n_initial_ppm >= 0.0 &&
        isfinite(cfg.volatiles.organic_n_initial_ppm)
    ) || throw(
        ArgumentError(
            "organic_n_initial_ppm must be >= 0 and finite, got $(cfg.volatiles.organic_n_initial_ppm)",
        ),
    )
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
    (cfg.volatiles.melt_feo_wtpct >= 0.0 && isfinite(cfg.volatiles.melt_feo_wtpct)) ||
        throw(
            ArgumentError(
                "melt_feo_wtpct must be >= 0 and finite, got $(cfg.volatiles.melt_feo_wtpct)",
            ),
        )
    (0.0 <= cfg.volatiles.x_sio2 <= 1.0 && isfinite(cfg.volatiles.x_sio2)) || throw(
        ArgumentError("x_sio2 must be in [0, 1] and finite, got $(cfg.volatiles.x_sio2)"),
    )
    (0.0 <= cfg.volatiles.x_al2o3 <= 1.0 && isfinite(cfg.volatiles.x_al2o3)) || throw(
        ArgumentError("x_al2o3 must be in [0, 1] and finite, got $(cfg.volatiles.x_al2o3)"),
    )
    (0.0 <= cfg.volatiles.x_tio2 <= 1.0 && isfinite(cfg.volatiles.x_tio2)) || throw(
        ArgumentError("x_tio2 must be in [0, 1] and finite, got $(cfg.volatiles.x_tio2)"),
    )

    # Escape checks
    (cfg.escape.M_planet > 0.0 && isfinite(cfg.escape.M_planet)) ||
        throw(ArgumentError("M_planet must be > 0 and finite, got $(cfg.escape.M_planet)"))
    (cfg.escape.R_planet > 0.0 && isfinite(cfg.escape.R_planet)) ||
        throw(ArgumentError("R_planet must be > 0 and finite, got $(cfg.escape.R_planet)"))
    if cfg.escape.active
        isapprox(cfg.escape.R_planet, cfg.geometry.rplanet; rtol=0.01) || throw(
            ArgumentError(
                "escape.R_planet ($(cfg.escape.R_planet)) must match simulated planet radius geometry.rplanet ($(cfg.geometry.rplanet)) when escape.active=true",
            ),
        )
    end
    (cfg.escape.T_exobase > 0.0 && isfinite(cfg.escape.T_exobase)) || throw(
        ArgumentError("T_exobase must be > 0 and finite, got $(cfg.escape.T_exobase)")
    )
    (cfg.escape.R_exobase >= cfg.escape.R_planet && isfinite(cfg.escape.R_exobase)) ||
        throw(
            ArgumentError(
                "R_exobase must be >= R_planet ($(cfg.escape.R_planet)) and finite, got $(cfg.escape.R_exobase)",
            ),
        )
    cfg.escape.species in Set([:H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2]) ||
        throw(
            ArgumentError(
                "escape species must be one of :H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2, got $(cfg.escape.species)",
            ),
        )

    # Core formation validation
    if cfg.coreformation.percolation_active || cfg.coreformation.settling_active
        cf = cfg.coreformation
        (isfinite(cf.rho_metal) && cf.rho_metal > cfg.materials.rhosolidm[1]) || throw(
            ArgumentError(
                "rho_metal ($(cf.rho_metal)) must be finite and exceed silicate rock density ($(cfg.materials.rhosolidm[1]))",
            ),
        )
        (cf.eta_metal > 0.0 && isfinite(cf.eta_metal)) || throw(
            ArgumentError("eta_metal must be positive and finite, got $(cf.eta_metal)")
        )
        (cf.k_metal > 0.0 && isfinite(cf.k_metal)) ||
            throw(ArgumentError("k_metal must be positive and finite, got $(cf.k_metal)"))
        (cf.rhocp_metal > 0.0 && isfinite(cf.rhocp_metal)) || throw(
            ArgumentError("rhocp_metal must be positive and finite, got $(cf.rhocp_metal)"),
        )
        (0.0 <= cf.Xfe_bulk <= 1.0) ||
            throw(ArgumentError("Xfe_bulk must be in [0, 1], got $(cf.Xfe_bulk)"))
        (0.0 < cf.phi_pack <= 1.0) ||
            throw(ArgumentError("phi_pack must be in (0, 1], got $(cf.phi_pack)"))
        (0.0 <= cf.phi_residual <= cf.phi_crit_perc < cf.phi_pack) || throw(
            ArgumentError(
                "phi_residual ($(cf.phi_residual)) must be <= phi_crit_perc ($(cf.phi_crit_perc)) < phi_pack ($(cf.phi_pack))",
            ),
        )
        (0.0 <= cf.F_settle_start <= cf.F_perc_end <= 1.0) || throw(
            ArgumentError(
                "F_settle_start ($(cf.F_settle_start)) must be <= F_perc_end ($(cf.F_perc_end)) in [0, 1]",
            ),
        )
        (0.0 < cf.cfl_settling <= 1.0) ||
            throw(ArgumentError("cfl_settling must be in (0, 1], got $(cf.cfl_settling)"))
        cf.max_subcycles > 0 ||
            throw(ArgumentError("max_subcycles must be > 0, got $(cf.max_subcycles)"))
        cf.droplet_size_mode in Set([:fixed, :weber_mean, :weber_turbulent]) || throw(
            ArgumentError(
                "droplet_size_mode must be one of :fixed, :weber_mean, :weber_turbulent, got $(cf.droplet_size_mode)",
            ),
        )
        (cf.droplet_diameter_fixed > 0.0 && isfinite(cf.droplet_diameter_fixed)) || throw(
            ArgumentError(
                "droplet_diameter_fixed must be positive, got $(cf.droplet_diameter_fixed)",
            ),
        )
        (cf.sigma_metal_silicate > 0.0 && isfinite(cf.sigma_metal_silicate)) || throw(
            ArgumentError(
                "sigma_metal_silicate must be positive, got $(cf.sigma_metal_silicate)"
            ),
        )
        (cf.We_crit > 0.0 && isfinite(cf.We_crit)) ||
            throw(ArgumentError("We_crit must be positive, got $(cf.We_crit)"))
        (cf.hindered_exponent >= 0.0 && isfinite(cf.hindered_exponent)) || throw(
            ArgumentError(
                "hindered_exponent must be non-negative, got $(cf.hindered_exponent)"
            ),
        )
    end

    if cfg.coreformation.percolation_active
        cf = cfg.coreformation
        cf.T_eutectic < cfg.melting.T_solidus[1] || throw(
            ArgumentError(
                "T_eutectic ($(cf.T_eutectic)) must be below silicate solidus ($(cfg.melting.T_solidus[1])) for percolation",
            ),
        )
        (cf.dT_metal > 0.0 && isfinite(cf.dT_metal)) ||
            throw(ArgumentError("dT_metal must be positive, got $(cf.dT_metal)"))
        (cf.k_metal_ref > 0.0 && isfinite(cf.k_metal_ref)) ||
            throw(ArgumentError("k_metal_ref must be positive, got $(cf.k_metal_ref)"))
        (cf.perm_exponent >= 0.0 && isfinite(cf.perm_exponent)) || throw(
            ArgumentError("perm_exponent must be non-negative, got $(cf.perm_exponent)")
        )
    end

    if cfg.coreformation.settling_active
        cfg.melting.active || throw(
            ArgumentError(
                "settling_active requires melting.active = true for silicate melt fraction",
            ),
        )
    end

    if cfg.volatiles.active
        @warn "VolatilesConfig active=true: multi-species H-C-N-S volatile solubility, gas speciation, and organic devolatilization operate as a standalone thermodynamic library; dynamic reactive transport is not yet coupled to the 2D Stokes-Darcy fluid flow solver."
    end

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
                kwargs[fname] = convert(ftype, val)
            elseif ftype === Symbol && val isa AbstractString
                kwargs[fname] = Symbol(val)
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
    "escape",
    "coreformation",
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
    esc = if haskey(parsed, "escape")
        parsed_esc = parsed["escape"]
        def_esc =
            if !haskey(parsed_esc, "R_planet") && geom.rplanet != def.geometry.rplanet
                EscapeConfig(
                    def.escape.active,
                    def.escape.M_planet,
                    geom.rplanet,
                    def.escape.T_exobase,
                    geom.rplanet,
                    def.escape.species,
                )
            else
                def.escape
            end
        _dict_to_struct(EscapeConfig, parsed_esc, def_esc)
    else
        if geom.rplanet != def.geometry.rplanet
            EscapeConfig(
                def.escape.active,
                def.escape.M_planet,
                geom.rplanet,
                def.escape.T_exobase,
                geom.rplanet,
                def.escape.species,
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
        escape=esc,
        coreformation=coreform,
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
        else
            d[String(fname)] = val
        end
    end
    return d
end

"""
Saves a `SimulationConfig` to an IO stream or a `.toml` file.

$(SIGNATURES)

# Arguments
- `io_or_path`: Output IO stream or file path.
- `cfg`: Configuration to serialize.
"""
function save_config(io::IO, cfg::SimulationConfig)
    d = Dict{String,Any}(
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
        "escape" => _struct_to_dict(cfg.escape),
        "coreformation" => _struct_to_dict(cfg.coreformation),
    )
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
