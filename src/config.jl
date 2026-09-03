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
    xsize::Float64 = 14_000.0
    ysize::Float64 = 14_000.0
    Nx::Int = 15
    Ny::Int = 15
end

"""
Planetary body geometry parameters.

$(FIELDS)
"""
Base.@kwdef struct GeometryConfig
    rplanet::Float64 = 5_000.0
    rcrust::Float64 = 4_800.0
    xcenter::Float64 = 7_000.0
    ycenter::Float64 = 7_000.0
    psurface::Float64 = 1.0e+3
end

"""
Timestepping and temporal integration parameters.

$(FIELDS)
"""
Base.@kwdef struct TimeConfig
    dt_initial::Float64 = 1.0e+11
    dt_longest::Float64 = 1.0e+11
    dtcoefdn::Float64 = 0.5
    dtcoefup::Float64 = 1.2
    dtstep::Int = 200
    dxymax::Float64 = 0.05
    vpratio::Float64 = 1.0 / 3.0
    DTmax::Float64 = 20.0
    yearlength::Float64 = 365.25 * 24 * 3600
    start_time::Float64 = 2.25e6 * (365.25 * 24 * 3600)
    endtime::Float64 = 15.0e6 * (365.25 * 24 * 3600)
    start_step::Int = 1
    n_steps::Int = 10
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
end

"""
Material phase properties for 3-phase staggered grid.

Index 1: Planetesimal core / mantle
Index 2: Porous silicate crust / rock
Index 3: Sticky air / space

$(FIELDS)
"""
Base.@kwdef struct MaterialConfig
    rhosolidm::SVector{3, Float64} = SVector{3, Float64}([3300.0, 3300.0, 1.0])
    rhofluidm::SVector{3, Float64} = SVector{3, Float64}([1000.0, 1000.0, 1.0])
    etasolidm::SVector{3, Float64} = SVector{3, Float64}([1.0e+19, 1.0e+19, 1.0e+16])
    etasolidmm::SVector{3, Float64} = SVector{3, Float64}([1.0e+19, 1.0e+19, 1.0e+16])
    etafluidm::SVector{3, Float64} = SVector{3, Float64}([1.0e+12, 1.0e+12, 1.0e-03])
    etafluidmm::SVector{3, Float64} = SVector{3, Float64}([1.0e-03, 1.0e-03, 1.0e-03])
    rhocpsolidm::SVector{3, Float64} = SVector{3, Float64}([3.3e+06, 3.3e+06, 3.0e+06])
    rhocpfluidm::SVector{3, Float64} = SVector{3, Float64}([1.0e+06, 1.0e+06, 3.0e+06])
    alphasolidm::SVector{3, Float64} = SVector{3, Float64}([3.0e-05, 3.0e-05, 0.0])
    alphafluidm::SVector{3, Float64} = SVector{3, Float64}([5.0e-05, 5.0e-05, 0.0])
    ksolidm::SVector{3, Float64} = SVector{3, Float64}([3.0, 3.0, 3000.0])
    kfluidm::SVector{3, Float64} = SVector{3, Float64}([50.0, 50.0, 3000.0])
    gggsolidm::SVector{3, Float64} = SVector{3, Float64}([1.0e+10, 1.0e+10, 1.0e+10])
    frictsolidm::SVector{3, Float64} = SVector{3, Float64}([0.6, 0.6, 0.0])
    cohessolidm::SVector{3, Float64} = SVector{3, Float64}([1.0e+08, 1.0e+08, 1.0e+08])
    tenssolidm::SVector{3, Float64} = SVector{3, Float64}([6.0e+07, 6.0e+07, 6.0e+07])
    kphim0::SVector{3, Float64} = SVector{3, Float64}([1.0e-13, 1.0e-13, 1.0e-17])
    tkm0::SVector{3, Float64} = SVector{3, Float64}([170.0, 170.0, 170.0])
end

"""
Output and storage parameters.

$(FIELDS)
"""
Base.@kwdef struct OutputConfig
    output_dir::String = "output"
    savematstep::Int = 10
    visstep::Int = 1
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
    materials::MaterialConfig = MaterialConfig()
    output::OutputConfig = OutputConfig()
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
    # Grid checks: must be valid dimensions and match compile-time stencils
    cfg.grid.Nx >= 3 || throw(ArgumentError("Grid Nx must be >= 3, got $(cfg.grid.Nx)"))
    cfg.grid.Ny >= 3 || throw(ArgumentError("Grid Ny must be >= 3, got $(cfg.grid.Ny)"))
    cfg.grid.xsize > 0.0 || throw(ArgumentError("Domain xsize must be > 0, got $(cfg.grid.xsize)"))
    cfg.grid.ysize > 0.0 || throw(ArgumentError("Domain ysize must be > 0, got $(cfg.grid.ysize)"))
    if cfg.grid.Nx != Nx || cfg.grid.Ny != Ny
        throw(ArgumentError("Configured grid size ($(cfg.grid.Nx)x$(cfg.grid.Ny)) does not match compiled grid size ($(Nx)x$(Ny)). Changing grid resolution requires updating constants and recompiling."))
    end
    if cfg.grid.xsize != xsize || cfg.grid.ysize != ysize
        throw(ArgumentError("Configured domain size ($(cfg.grid.xsize)x$(cfg.grid.ysize)) does not match compiled domain size ($(xsize)x$(ysize))."))
    end

    # Geometry checks
    cfg.geometry.rplanet > 0.0 || throw(ArgumentError("Planet radius must be > 0, got $(cfg.geometry.rplanet)"))
    cfg.geometry.rcrust > 0.0 || throw(ArgumentError("Crust radius must be > 0, got $(cfg.geometry.rcrust)"))
    cfg.geometry.rcrust <= cfg.geometry.rplanet || throw(ArgumentError("Crust radius must be <= planet radius"))
    if cfg.geometry.rplanet != rplanet || cfg.geometry.rcrust != rcrust
        throw(ArgumentError("Configured geometry (rplanet=$(cfg.geometry.rplanet), rcrust=$(cfg.geometry.rcrust)) does not match compiled geometry (rplanet=$rplanet, rcrust=$rcrust)."))
    end

    # Time checks
    cfg.time.dt_initial > 0.0 || throw(ArgumentError("Initial dt must be > 0, got $(cfg.time.dt_initial)"))
    cfg.time.dt_longest >= cfg.time.dt_initial || throw(ArgumentError("dt_longest must be >= dt_initial"))
    cfg.time.n_steps >= 1 || throw(ArgumentError("n_steps must be >= 1, got $(cfg.time.n_steps)"))
    cfg.time.start_step >= 1 || throw(ArgumentError("start_step must be >= 1, got $(cfg.time.start_step)"))
    isfinite(cfg.time.dt_initial) || throw(ArgumentError("dt_initial must be finite"))
    isfinite(cfg.time.dt_longest) || throw(ArgumentError("dt_longest must be finite"))

    # Poroelasticity checks
    cfg.poroelasticity.betasolid >= 0.0 || throw(ArgumentError("betasolid must be >= 0, got $(cfg.poroelasticity.betasolid)"))
    cfg.poroelasticity.betafluid >= 0.0 || throw(ArgumentError("betafluid must be >= 0, got $(cfg.poroelasticity.betafluid)"))
    isfinite(cfg.poroelasticity.betasolid) || throw(ArgumentError("betasolid must be finite"))
    isfinite(cfg.poroelasticity.betafluid) || throw(ArgumentError("betafluid must be finite"))
    0.0 < cfg.poroelasticity.phimin < cfg.poroelasticity.phimax < 1.0 || throw(
        ArgumentError("Porosity bounds must satisfy 0 < phimin < phimax < 1, got phimin=$(cfg.poroelasticity.phimin), phimax=$(cfg.poroelasticity.phimax)")
    )

    # Solver checks
    cfg.solver.titermax >= 1 || throw(ArgumentError("titermax must be >= 1, got $(cfg.solver.titermax)"))
    cfg.solver.nplast >= 1 || throw(ArgumentError("nplast must be >= 1, got $(cfg.solver.nplast)"))
    cfg.solver.titermax <= cfg.solver.nplast || throw(ArgumentError("titermax ($(cfg.solver.titermax)) must be <= nplast ($(cfg.solver.nplast)) to prevent array bounds overflow in plastic convergence tracking"))
    cfg.solver.etamin > 0.0 || throw(ArgumentError("etamin must be > 0, got $(cfg.solver.etamin)"))
    cfg.solver.etamax >= cfg.solver.etamin || throw(ArgumentError("etamax must be >= etamin"))
    cfg.solver.etaphikoef > 0.0 || throw(ArgumentError("etaphikoef must be > 0, got $(cfg.solver.etaphikoef)"))

    # Output checks
    cfg.output.savematstep >= 1 || throw(ArgumentError("savematstep must be >= 1, got $(cfg.output.savematstep)"))
    cfg.output.visstep >= 1 || throw(ArgumentError("visstep must be >= 1, got $(cfg.output.visstep)"))

    # Thermodynamics checks
    0.0 <= cfg.thermodynamics.ratio_al <= 1.0 || throw(ArgumentError("ratio_al must be in [0, 1], got $(cfg.thermodynamics.ratio_al)"))
    0.0 <= cfg.thermodynamics.ratio_fe <= 1.0 || throw(ArgumentError("ratio_fe must be in [0, 1], got $(cfg.thermodynamics.ratio_fe)"))
    cfg.thermodynamics.tmfluidphase < cfg.thermodynamics.tmsolidphase || throw(ArgumentError("tmfluidphase ($(cfg.thermodynamics.tmfluidphase)) must be < tmsolidphase ($(cfg.thermodynamics.tmsolidphase))"))
    cfg.thermodynamics.Lᶠ > 0.0 || throw(ArgumentError("Lᶠ must be > 0, got $(cfg.thermodynamics.Lᶠ)"))
    cfg.thermodynamics.E_al > 0.0 && isfinite(cfg.thermodynamics.E_al) || throw(ArgumentError("E_al must be > 0 and finite"))
    cfg.thermodynamics.f_al > 0.0 && isfinite(cfg.thermodynamics.f_al) || throw(ArgumentError("f_al must be > 0 and finite"))
    cfg.thermodynamics.t_half_al > 0.0 && isfinite(cfg.thermodynamics.t_half_al) || throw(ArgumentError("t_half_al must be > 0 and finite"))
    cfg.thermodynamics.E_fe > 0.0 && isfinite(cfg.thermodynamics.E_fe) || throw(ArgumentError("E_fe must be > 0 and finite"))
    cfg.thermodynamics.f_fe > 0.0 && isfinite(cfg.thermodynamics.f_fe) || throw(ArgumentError("f_fe must be > 0 and finite"))
    cfg.thermodynamics.t_half_fe > 0.0 && isfinite(cfg.thermodynamics.t_half_fe) || throw(ArgumentError("t_half_fe must be > 0 and finite"))

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

    if cfg.materials.rhosolidm != rhosolidm || cfg.materials.rhofluidm != rhofluidm ||
       cfg.materials.etasolidm != etasolidm || cfg.materials.etasolidmm != etasolidmm ||
       cfg.materials.etafluidm != etafluidm || cfg.materials.etafluidmm != etafluidmm ||
       cfg.materials.ksolidm != ksolidm || cfg.materials.kfluidm != kfluidm
        throw(ArgumentError("Overriding [materials] arrays at runtime is not supported because marker property assignment uses compiled constants. Recompilation is required to modify material properties."))
    end

    return nothing
end

"""
Helper function to convert TOML-parsed dictionary into a typed struct with defaults.
"""
function _dict_to_struct(::Type{T}, d::Dict{String, Any}, defaults::T) where {T}
    # Check for unknown / misspelled keys
    for k in keys(d)
        if !hasfield(T, Symbol(k))
            throw(ArgumentError("Unknown configuration key '$k' in [$(nameof(T))] section."))
        end
    end

    kwargs = Dict{Symbol, Any}()
    for fname in fieldnames(T)
        sname = String(fname)
        ftype = fieldtype(T, fname)
        if haskey(d, sname)
            val = d[sname]
            if ftype <: SVector
                expected_len = length(ftype)
                if !(val isa AbstractVector) || length(val) != expected_len
                    throw(ArgumentError("Field '$sname' in [$(nameof(T))] must have exactly $expected_len elements, got $(val)"))
                end
                kwargs[fname] = SVector{expected_len, eltype(ftype)}(val)
            elseif ftype <: Real && !(val isa ftype)
                kwargs[fname] = convert(ftype, val)
            else
                kwargs[fname] = val
            end
        else
            kwargs[fname] = getfield(defaults, fname)
        end
    end
    return T(; kwargs...)
end

const VALID_SECTIONS = Set(["grid", "geometry", "time", "solver", "poroelasticity", "thermodynamics", "materials", "output"])

"""
Loads, merges, and validates a `SimulationConfig` from a TOML file or string.

$(SIGNATURES)

# Arguments
- `source`: File path to `.toml` file, or raw TOML string.

# Returns
- `cfg::SimulationConfig`: Validated configuration struct.
"""
function load_config(source::AbstractString)::SimulationConfig
    # Resolve file path: direct path or relative to Erebus package root
    resolved_path = if isfile(source)
        source
    elseif isfile(joinpath(@__DIR__, "..", source))
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

    grid = haskey(parsed, "grid") ? _dict_to_struct(GridConfig, parsed["grid"], def.grid) : def.grid
    geom = haskey(parsed, "geometry") ? _dict_to_struct(GeometryConfig, parsed["geometry"], def.geometry) : def.geometry
    time = haskey(parsed, "time") ? _dict_to_struct(TimeConfig, parsed["time"], def.time) : def.time
    solv = haskey(parsed, "solver") ? _dict_to_struct(SolverConfig, parsed["solver"], def.solver) : def.solver
    poro = haskey(parsed, "poroelasticity") ? _dict_to_struct(PoroelasticConfig, parsed["poroelasticity"], def.poroelasticity) : def.poroelasticity
    therm = haskey(parsed, "thermodynamics") ? _dict_to_struct(ThermalConfig, parsed["thermodynamics"], def.thermodynamics) : def.thermodynamics
    mat = haskey(parsed, "materials") ? _dict_to_struct(MaterialConfig, parsed["materials"], def.materials) : def.materials
    out = haskey(parsed, "output") ? _dict_to_struct(OutputConfig, parsed["output"], def.output) : def.output

    cfg = SimulationConfig(
        grid = grid,
        geometry = geom,
        time = time,
        solver = solv,
        poroelasticity = poro,
        thermodynamics = therm,
        materials = mat,
        output = out
    )

    validate_config(cfg)
    return cfg
end

"""
Converts a struct into a Dict suitable for TOML serialization.
"""
function _struct_to_dict(s)
    d = Dict{String, Any}()
    for fname in fieldnames(typeof(s))
        val = getfield(s, fname)
        if val isa SVector
            d[String(fname)] = collect(val)
        else
            d[String(fname)] = val
        end
    end
    return d
end

"""
Saves a `SimulationConfig` to a `.toml` file.

$(SIGNATURES)

# Arguments
- `path`: Output file path.
- `cfg`: Configuration to serialize.
"""
function save_config(path::AbstractString, cfg::SimulationConfig)
    d = Dict{String, Any}(
        "grid" => _struct_to_dict(cfg.grid),
        "geometry" => _struct_to_dict(cfg.geometry),
        "time" => _struct_to_dict(cfg.time),
        "solver" => _struct_to_dict(cfg.solver),
        "poroelasticity" => _struct_to_dict(cfg.poroelasticity),
        "thermodynamics" => _struct_to_dict(cfg.thermodynamics),
        "materials" => _struct_to_dict(cfg.materials),
        "output" => _struct_to_dict(cfg.output),
    )
    open(path, "w") do io
        TOML.print(io, d; sorted=true)
    end
    return path
end
