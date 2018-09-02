# Unit test for the lin-solve mehtods (in src/LinSolver.jl)
# The "Gun" problem form gun_native.jl

# Intended to be run from nep-pack/ directory or nep-pack/test directory
push!(LOAD_PATH, string(@__DIR__, "/../src"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using LinearAlgebra
using Random
using Test

@testset "linsolvers" begin
    TOL = 1e-10;
    nep = nep_gallery("nlevp_native_gun")

    n = size(nep,1);
    λ = 250^2+1im;

    x=ones(n);

    linsolver1=backslash_linsolvercreator(nep,λ);
    y1=lin_solve(linsolver1,x)

    linsolver2=default_linsolvercreator(nep,λ);
    y2=lin_solve(linsolver1,x)

    @test y1 ≈ y2

    y3=compute_Mder(nep,λ)\x;
    @test y1 ≈ y3
    @test y2 ≈ y3


    println("norm(y1,1)=",norm(y1));
    r=compute_Mder(nep,λ)*y1-x;
    println("norm(residual)=",norm(r));
    println("x'*r=",x'*r);

end
