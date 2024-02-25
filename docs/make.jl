using Documenter
using VectorEnumeration

makedocs(
    sitename = "VectorEnumeration.jl",
    authors = "Florian Heine and Nicolaus Jacobsen",
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
    repo = "github.com/Ktrompfl/VectorEnumeration.jl.git"
)
