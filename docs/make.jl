using Pkg
Pkg.develop(PackageSpec(; path=joinpath(@__DIR__, "..")))

using Documenter
using Erebus

DocMeta.setdocmeta!(Erebus, :DocTestSetup, :(using Erebus); recursive=true)

makedocs(;
    modules=[Erebus],
    authors="Tim Lichtenberg and Forming Worlds Lab contributors",
    sitename="Erebus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://formingworlds.github.io/Erebus.jl/stable/",
        edit_link="main",
        size_threshold_warn=800 * 1024,
        size_threshold=1200 * 1024,
        search_size_threshold_warn=1200 * 1024,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Quickstart" => "tutorials/quickstart.md",
            "2D Hydrothermal Circulation" => "tutorials/hydrothermal_circulation.md",
            "Planetesimal Differentiation" => "tutorials/planetesimal_differentiation.md",
            "Growth to Lunar Mass" => "tutorials/growth_to_lunar_mass.md",
        ],
        "How-To Guides" => [
            "Installation" => "howto/installation.md",
            "Configuration (.toml)" => "howto/configuration.md",
            "Running Simulations" => "howto/running.md",
            "Distributed Simulations (MPI)" => "howto/distributed_simulations.md",
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
            "Planetesimal Accretion Mechanics" => "explanations/accretion_mechanics.md",
            "Telescoping Domain Dynamics" => "explanations/telescoping_domain.md",
            "Magma Ascent & Matrix Compaction" => "explanations/magma_compaction.md",
        ],
        "Reference" => [
            "Configuration Schema" => "reference/config_schema.md",
            "API Reference" => "reference/api.md",
            "Physical Benchmarks" => [
                "Overview" => "reference/benchmarks/index.md",
                "1D Terzaghi Consolidation" => "reference/benchmarks/terzaghi_consolidation.md",
                "Thermal Conduction & Geometry" => "reference/benchmarks/thermal_conduction.md",
                "Permeability & Hydrofracture" => "reference/benchmarks/permeability.md",
                "Radionuclide Decay" => "reference/benchmarks/radionuclides.md",
                "Fluid Viscosity & Phase Changes" => "reference/benchmarks/fluid_viscosity.md",
                "Surface Radiation & Disk Evolution" => "reference/benchmarks/disk_radiation.md",
                "Hydrothermal Reactions" => "reference/benchmarks/hydrothermal_reactions.md",
                "Silicate Rock Melting" => "reference/benchmarks/rock_melting.md",
                "Cold Surface Venting" => "reference/benchmarks/cold_surface_venting.md",
                "Cold Lid Hydrofracture & Ice Sealing" => "reference/benchmarks/hydrofracture_venting.md",
                "Dehydration-Darcy Coupling & Venting" => "reference/benchmarks/dehydration_darcy_coupling.md",
                "H-C-N-S Volatile Solubility & Speciation" => "reference/benchmarks/hcns_solubility.md",
                "Volatile Retention Floors & Vent Drainage" => "reference/benchmarks/volatile_retention.md",
                "Redox Buffers & Electron Accounting" => "reference/benchmarks/redox_buffers.md",
                "Jeans Kinetic Atmospheric Escape" => "reference/benchmarks/jeans_escape.md",
                "Iron Core Formation" => "reference/benchmarks/core_formation.md",
                "Core Geochemistry & Volatiles" => "reference/benchmarks/core_geochemistry.md",
                "Normative Accessory Minerals" => "reference/benchmarks/mineral_assemblages.md",
                "Hydrothermal Subgrid Convection" => "reference/benchmarks/hydrothermal_convection.md",
                "Planetesimal Accretion & Impact Heating" => "reference/benchmarks/planetesimal_accretion.md",
                "Telescoping Computational Domain" => "reference/benchmarks/telescoping_domain.md",
                "Volatile & Refractory Mixtures" => "reference/benchmarks/volatile_mixtures.md",
                "Multi-Stage Accretion Sequence" => "reference/benchmarks/multistage_accretion.md",
                "Coupled Atmosphere & Gas Envelopes" => "reference/benchmarks/coupled_atmosphere.md",
                "Magma Compaction & Decompression Exsolution" => "reference/benchmarks/magma_compaction.md",
                "Crustal Sill Cooling & Magma-Hydrothermal Coupling" => "reference/benchmarks/sill_cooling.md",
            ],
            "Bibliography" => "reference/bibliography.md",
        ],
        "Community" => [
            "Contributing" => "community/contributing.md",
            "Acknowledgements" => "community/acknowledgements.md",
        ],
    ],
    warnonly=[:cross_references],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(;
        repo="github.com/FormingWorlds/Erebus.jl.git",
        devbranch="main",
        push_preview=true,
        versions=["stable" => "dev", "dev" => "dev"],
    )
end
