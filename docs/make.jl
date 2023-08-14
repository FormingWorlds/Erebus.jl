push!(LOAD_PATH, "../src/")

using Documenter
using Erebus

makedocs(
    sitename = "Erebus.jl",
    modules = [Erebus],
    authors="Beat Hubmann",
    repo="https://github.com/FormingWorlds/Erebus.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FormingWorlds/Erebus.jl",
    devbranch="main",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
