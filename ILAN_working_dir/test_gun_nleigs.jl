using NonlinearEigenproblems, Random, SparseArrays, PyPlot, LinearAlgebra

println("Loading the:", `gun`,"problem")
nep = nep_gallery("nlevp_native_gun")
println("Shifting and scaling the problem")
include("shift_and_scale_gun.jl");


# run tiar
println("Run nleigs")
θ=range(0,stop=π,length=100); r=50000;
Σ=r*cos.(θ) + 1im*r*sin.(θ)

# define the set of pole candidates
Ξ = -10 .^ range(-8, stop = 8, length = 10000) .- r^2

λ=nleigs(nep,Σ;Ξ=Ξ,logger=2)
