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

TIAR=@testset "TIAR" begin
    @testset "accuracy eigenpairs" begin
        (λ,Q)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
        @testset "TIAR eigval[$i]" for i in 1:length(λ)
            @test norm(compute_Mlincomb(dep,λ[i],Q[:,i]))<eps()*100;
        end
    end

    @testset "orthogonalization" begin

    # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS
    @testset "DGKS" begin
        (λ,Q,err,Z)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
        @test norm(Z'*Z-eye(size(Z,2)))<1e-6
     end

     @testset "User provided doubleGS" begin
         (λ,Q,err,Z)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
         @test norm(Z'*Z-eye(size(Z,2)))<1e-6
      end

      @testset "ModifiedGramSchmidt" begin
          (λ,Q,err,Z)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
          @test norm(Z'*Z-eye(size(Z,2)))<1e-6
      end

       @testset "ClassicalGramSchmidt" begin
           (λ,Q,err,Z)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
           @test norm(Z'*Z-eye(size(Z,2)))<1e-6
       end
    end

    # iar and tiar ara mathematically equivalent it maxit<<nep
    # verify the equivalence numerically
    @testset "equivalance with IAR" begin
        (λ_tiar,Q_tiar)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
        (λ_iar,Q_iar)=tiar(dep,σ=2.0,γ=3,Neig=4,v=ones(n),displaylevel=0,maxit=50,tol=eps()*100);
        @test norm(λ_tiar-λ_iar)<1e-6
        @test norm(Q_tiar-Q_iar)<1e-6
    end


end
