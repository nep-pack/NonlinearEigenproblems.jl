push!(LOAD_PATH, string(@__DIR__, "/../src"))
using Documenter, NEPCore, NEPTypes, NEPSolver


makedocs(
    clean = true,
    doctest = false,
    pages = Any[
        "Home" => "index.md",
        "NEP Methods" => "methods.md",
        "NEP Transformations" => "transformations.md",
        "NEP Gallery" => "gallery.md"                
    ]
)


