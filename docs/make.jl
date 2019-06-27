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
        "NEP Types" => "types.md",
        "LinSolver" => "linsolvers.md",
        "Error measure" => "errmeasure.md",
        "NEP Transformations" => "transformations.md",
        "NEP Gallery" => "gallery.md",
        "Tutorial 1 (ABC)" => "movebc_tutorial.md",
        "Tutorial 2 (BEM)" => "bemtutorial.md",
        "Tutorial 3 (Deflation)" => "deflate_tutorial.md",
        "Tutorial 4 (Python NEP)" => "tutorial_call_python.md",
        "Tutorial 5 (MATLAB 1)" => "tutorial_matlab1.md",
        "Tutorial 6 (FORTRAN 1)" => "tutorial_fortran1.md",
        "Tutorial 7 (Hydrodynamic stability PEP)" => "hydrotutorial.md"
    ]
)


#deploydocs(
#    repo = "github.com/nep-pack/NonlinearEigenproblems.jl.git",
#    target = "build",
#)
