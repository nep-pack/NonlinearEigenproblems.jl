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
        "Logger" => "logger.md",
        "NEP Transformations" => "transformations.md",
        "NEP Gallery" => "gallery.md",
        "Tutorial 1 (ABC)" => "movebc_tutorial.md",
        "Tutorial 2 (Contour)" => "tutorial_contour.md",
        "Tutorial 3 (BEM)" => "bemtutorial.md",
        "Tutorial 4 (Deflation)" => "deflate_tutorial.md",
        "Tutorial 5 (Python NEP)" => "tutorial_call_python.md",
        "Tutorial 6 (MATLAB 1)" => "tutorial_matlab1.md",
        "Tutorial 7 (FORTRAN 1)" => "tutorial_fortran1.md",
        "Tutorial 8 (gmsh + nanophotonics)" => "tutorial_nano1.md",
        "Tutorial 9 (New solver)" => "tutorial_newmethod.md",
        "Tutorial 10 (Orr-Somerfeld)" => "hydrotutorial.md",
    ]
)


#deploydocs(
#    repo = "github.com/nep-pack/NonlinearEigenproblems.jl.git",
#    target = "build",
#)
