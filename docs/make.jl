using Pkg
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))

using Documenter
using Erebus

DocMeta.setdocmeta!(Erebus, :DocTestSetup, :(using Erebus); recursive=true)

makedocs(
    modules = [Erebus, Erebus.Config, Erebus.Geometry, Erebus.Physics, Erebus.Particles, Erebus.Numerics, Erebus.Simulation],
    authors = "Tim Lichtenberg and Forming Worlds Lab contributors",
    sitename = "Erebus.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://formingworlds.github.io/Erebus.jl/stable/",
        edit_link = "main",
        size_threshold_warn = 250 * 1024,
        size_threshold = 350 * 1024,
    ),
    pages = [
        "Home" => "index.md",

        "Tutorials" => [
            "Quickstart" => "tutorials/quickstart.md",
            "1D Terzaghi Consolidation" => "tutorials/terzaghi_consolidation.md",
            "2D Hydrothermal Circulation" => "tutorials/hydrothermal_circulation.md",
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
            "Verification & Benchmarks" => "explanations/verification.md",
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
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/FormingWorlds/Erebus.jl.git",
    devbranch = "main",
    push_preview = true,
)
