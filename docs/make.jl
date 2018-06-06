#push!(LOAD_PATH, string(@__DIR__, "/../src"))
using Documenter;
using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.NEPTypes
using NonlinearEigenproblems.NEPCore
#using NonlinearEigenproblems: NEPCore, NEPTypes, NEPSolver
#importall  NonlinearEigenproblems: NEPCore, NEPTypes, NEPSolver
#importall NonlinearEigenproblems

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


