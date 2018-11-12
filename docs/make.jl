push!(LOAD_PATH, string(@__DIR__, "/.."))
using Documenter;
using NonlinearEigenproblems

makedocs(
    clean = true,
    doctest = false,
    sitename = "NEP-PACK",
    pages = Any[
        "Home" => "index.md",
        "NEP Methods" => "methods.md",
        "LinSolver" => "linsolvers.md",
        "NEP Transformations" => "transformations.md",
        "NEP Gallery" => "gallery.md"
    ]
)


#deploydocs(
#    repo = "github.com/nep-pack/NonlinearEigenproblems.jl.git",
#    target = "build",
#)
