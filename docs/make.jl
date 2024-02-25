using Documenter
using VectorEnumeration

makedocs(
    sitename = "VectorEnumeration.jl",
    authors = "Florian Heine <fheine@rptu.de> and Nicolaus Jacobsen <jacobsen@rptu.de>",
    # modules = [VectorEnumeration],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    doctest = true,
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    warnonly = [:missing_docs]
)

deploydocs(
    repo = "https://github.com/Ktrompfl/VectorEnumeration.jl.git"
)
