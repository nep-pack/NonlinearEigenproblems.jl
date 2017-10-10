push!(LOAD_PATH, string(@__DIR__, "/../src"))
using Documenter, NEPCore, NEPTypes
#deploydocs(
#    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
#)

makedocs(
    clean = true,
    doctest = false,
)


