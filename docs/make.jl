using Pkg
Pkg.develop(PackageSpec(; path=joinpath(@__DIR__, "..")))

using Documenter
using Erebus

DocMeta.setdocmeta!(Erebus, :DocTestSetup, :(using Erebus); recursive=true)

makedocs(;
    modules=[
        Erebus,
        Erebus.Config,
        Erebus.Geometry,
        Erebus.Physics,
        Erebus.Particles,
        Erebus.Numerics,
        Erebus.Simulation,
    ],
    authors="Tim Lichtenberg and Forming Worlds Lab contributors",
    sitename="Erebus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://proteus-framework.org/Erebus.jl/stable/",
        edit_link="main",
        size_threshold_warn=350 * 1024,
        size_threshold=450 * 1024,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Quickstart" => "tutorials/quickstart.md",
            "1D Terzaghi Consolidation" => "tutorials/terzaghi_consolidation.md",
            "2D Hydrothermal Circulation" => "tutorials/hydrothermal_circulation.md",
            "Planetesimal Differentiation" => "tutorials/planetesimal_differentiation.md",
        ],
        "How-To Guides" => [
            "Installation" => "howto/installation.md",
            "Configuration (.toml)" => "howto/configuration.md",
            "Running Simulations" => "howto/running.md",
            "Parameter Exploration" => "howto/parameter_exploration.md",
            "Outputs & Checkpoints" => "howto/outputs.md",
        ],
        "Explanations" => [
            "Model Overview" => "explanations/model_overview.md",
            "Governing Equations" => "explanations/governing_equations.md",
            "Discretization & Numerics" => "explanations/discretization_numerics.md",
            "Silicate Melting & Soft Turbulence" => "explanations/rock_melting.md",
            "Protoplanetary Disk Evolution" => "explanations/disk_temperature_evolution.md",
            "Degassing & Cold Venting" => "explanations/degassing_and_venting.md",
            "Iron Core Formation & Metal Segregation" => "explanations/core_formation.md",
            "Verification & Benchmarks" => "explanations/verification.md",
        ],
        "Validation" => [
            "Overview" => "validation/index.md",
            "Thermal Conduction & Geometry" => "validation/thermal_conduction.md",
            "Permeability & Hydrofracture" => "validation/permeability.md",
            "Radionuclide Decay" => "validation/radionuclides.md",
            "Fluid Viscosity & Phase Changes" => "validation/fluid_viscosity.md",
            "Surface Radiation & Disk Evolution" => "validation/disk_radiation.md",
            "Hydrothermal Reactions" => "validation/hydrothermal_reactions.md",
            "Silicate Rock Melting" => "validation/rock_melting.md",
            "Cold Surface Venting" => "validation/cold_surface_venting.md",
            "Cold Lid Hydrofracture & Ice Sealing" => "validation/hydrofracture_venting.md",
            "H-C-N-S Volatile Solubility & Speciation" => "validation/hcns_solubility.md",
            "Volatile Retention Floors & Vent Drainage" => "validation/volatile_retention.md",
            "Jeans Kinetic Atmospheric Escape" => "validation/jeans_escape.md",
            "Iron Core Formation" => "validation/core_formation.md",
        ],
        "Reference" => [
            "Configuration Schema" => "reference/config_schema.md",
            "API Reference" => "reference/api.md",
            "Bibliography" => "reference/bibliography.md",
        ],
        "Community" => [
            "Contributing" => "community/contributing.md",
            "Acknowledgements" => "community/acknowledgements.md",
        ],
    ],
    warnonly=[:cross_references],
)

deploydocs(;
    repo="github.com/FormingWorlds/Erebus.jl.git",
    devbranch="main",
    push_preview=true,
    versions=["stable" => "dev", "dev" => "dev"],
)
