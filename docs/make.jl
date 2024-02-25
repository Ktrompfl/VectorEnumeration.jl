using Documenter
using VectorEnumeration

makedocs(
    sitename = "VectorEnumeration.jl",
    authors = "Florian Heine <fheine@rptu.de> and Nicolaus Jacobsen <jacobsen@rptu.de>",
    # modules = [VectorEnumeration],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    remotes = nothing,
    doctest = true,
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    warnonly = [:missing_docs]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
