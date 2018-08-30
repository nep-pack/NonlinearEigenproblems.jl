#  A Polynomial eigenvalue problem
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinearAlgebra
using Random

n=5;
Random.seed!(0) # reset the random seed
A0=randn(n,n);
A1=randn(n,n);
A2=randn(n,n);

println("Create a REP")
nep=REP([A0,A1,A2],[0,1,1im])

v=ones(n,1);
λ=1;
println("A test call to compute_Mlincomb")
z=compute_Mlincomb(nep,λ,v);
println("Running augnewton")
λ,x =augnewton(nep,displaylevel=1);
println("Computed solution to resnorm:",norm(compute_Mlincomb(nep,λ,x)))
