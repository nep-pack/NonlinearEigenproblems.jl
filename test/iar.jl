#Intended to be run from nep-pack/ directory or nep-pack/test directory
push!(LOAD_PATH, string(@__DIR__, "/../src"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using IterativeSolvers
using LinearAlgebra
using Test

#import NEPSolver.iar;
#include("../src/method_iar.jl");


# The user can create his own orthogonalization function to use in IAR
function doubleGS_function!(VV, vv, h)
    h[:]=VV'*vv; vv[:]=vv-VV*h; g=VV'*vv; vv[:]=vv-VV*g;
    h[:] = h[:]+g[:]; β=norm(vv); vv[:]=vv/β; return β
end
# Then it is needed to create a type to access to this function
abstract type DoubleGS <: IterativeSolvers.OrthogonalizationMethod end
# And then introduce a function dispatch for this new type in order to use
# the defined orthogonalization function
import IterativeSolvers.orthogonalize_and_normalize!
function orthogonalize_and_normalize!(V,v,h,::Type{DoubleGS})
    doubleGS_function!(V, v, h) end

@testset "IAR" begin
    dep=nep_gallery("dep0");
    n=size(dep,1);

    @testset "accuracy eigenpairs" begin
        (λ,Q)=iar(dep,σ=3,Neig=5,v=ones(n),
                  displaylevel=0,maxit=100,tol=eps()*100);
        @testset "IAR eigval[$i]" for i in 1:length(λ)
            @test norm(compute_Mlincomb(dep,λ[i],Q[:,i]))<eps()*100;
        end
    end

    @testset "orthogonalization" begin

    # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS
    @testset "DGKS" begin
        (λ,Q,err,V)=iar(dep,orthmethod=DGKS,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100)
        @test opnorm(V'*V - I) < 1e-6
     end

     @testset "User provided doubleGS" begin
         (λ,Q,err,V)=iar(dep,orthmethod=DoubleGS,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100)
         @test opnorm(V'*V - I) < 1e-6
      end

      @testset "ModifiedGramSchmidt" begin
          (λ,Q,err,V)=iar(dep,orthmethod=ModifiedGramSchmidt,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100)
          @test opnorm(V'*V - I) < 1e-6
      end

       @testset "ClassicalGramSchmidt" begin
           (λ,Q,err,V)=iar(dep,orthmethod=ClassicalGramSchmidt,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100)
           @test opnorm(V'*V - I) < 1e-6
       end
    end

end
