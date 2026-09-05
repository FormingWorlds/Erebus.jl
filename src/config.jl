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
