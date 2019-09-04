using NonlinearEigenproblems, Random, LinearAlgebra, PyPlot, Revise
include("extract_nleigs_structure.jl")


nep=nep_gallery("dep0");
unit_square = float([1+1im, 1-1im, -1-1im,-1+1im])
D, β, ξ, σ=nleigs_structure(nep,unit_square,);
