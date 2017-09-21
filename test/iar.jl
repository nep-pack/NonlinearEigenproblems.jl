# Run tests for the dep_distributed example

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/gallery_extra")
push!(LOAD_PATH, pwd()*"/src/gallery_extra/waveguide")

using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers
using IterativeSolvers


using Base.Test



dep=nep_gallery("dep0");
n=size(dep,1);

IAR=@testset "IAR" begin
    @testset "accuracy eigenpairs" begin

        (λ,Q)=iar(dep,orthmethod=ModifiedGramSchmidt,σ=3,Neig=5,v=ones(n),
            displaylevel=0,maxit=100,tol=eps()*100)

        @testset "IAR eigval[$i]" for i in 1:length(λ)
            @test norm(compute_Mlincomb(dep,λ[i],Q[:,i]))<eps()*100
        end
    end


end
