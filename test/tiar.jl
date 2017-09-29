# Run tests for the dep_distributed example

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))
using NEPSolver

using NEPCore
using NEPTypes
using Gallery
using LinSolvers
using IterativeSolvers
using Base.Test



# The user can create his own orthogonalization function to use in IAR
function doubleGS_function!(VV, vv, h)
    h[:]=VV'*vv; vv[:]=vv-VV*h; g=VV'*vv; vv[:]=vv-VV*g;
    h[] = h[]+g[]; β=norm(vv); vv[:]=vv/β; return β
end
# Then it is needed to create a type to access to this function
abstract type DoubleGS <: IterativeSolvers.OrthogonalizationMethod end
# And then introduce a function dispatch for this new type in order to use
# the defined orthogonalization function
import IterativeSolvers.orthogonalize_and_normalize!
function orthogonalize_and_normalize!(V,v,h,::Type{DoubleGS})
    doubleGS_function!(V, v, h) end

n=100;
dep=nep_gallery("dep0",n);
n=size(dep,1);

TIAR=@testset "TIAR" begin
    @testset "accuracy eigenpairs" begin
        (λ,Q)=tiar(dep,σ=0,Neig=3,v=ones(n),displaylevel=1,maxit=50,tol=eps()*100);
        @testset "TIAR eigval[$i]" for i in 1:length(λ)
            @test norm(compute_Mlincomb(dep,λ[i],Q[:,i]))<eps()*100;
        end
    end


end
