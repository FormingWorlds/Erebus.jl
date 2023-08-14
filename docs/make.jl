push!(LOAD_PATH, "../src/")

using Documenter
using Erebus

makedocs(
    sitename = "Erebus",
    format = Documenter.HTML(),
    modules = [Erebus]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
