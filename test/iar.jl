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

using Base.Test



dep=nep_gallery("dep0");
n=size(dep,1);

@testset "IAR" begin

    (λ,V)=iar(dep,σ=3,Neig=5,v=ones(n),
          displaylevel=1,maxit=100,tolerance=eps()*100)

    @testset "IAR eigval[$i]" for i in 1:length(λ)
      @test norm(compute_Mlincomb(dep,λ[i],V[:,i]))<eps()*100
    end

end
