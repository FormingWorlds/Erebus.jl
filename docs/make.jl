push!(LOAD_PATH, "../src/")

using Documenter
using Erebus

makedocs(
    modules = [Erebus],
    sitename = "Erebus.jl",
)

deploydocs(
    repo="github.com/FormingWorlds/Erebus.jl",
)
