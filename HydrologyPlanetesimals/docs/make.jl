using HydrologyPlanetesimals
using Documenter

DocMeta.setdocmeta!(HydrologyPlanetesimals, :DocTestSetup, :(using HydrologyPlanetesimals); recursive=true)

makedocs(;
    modules=[HydrologyPlanetesimals],
    authors="Beat Hubmann",
    repo="https://github.com/BeatHubmann/HydrologyPlanetesimals.jl/blob/{commit}{path}#{line}",
    sitename="HydrologyPlanetesimals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://BeatHubmann.github.io/HydrologyPlanetesimals.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/BeatHubmann/HydrologyPlanetesimals.jl",
    devbranch="main",
)
