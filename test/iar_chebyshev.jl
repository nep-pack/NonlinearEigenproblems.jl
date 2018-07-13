#Run tests for the dep_distributed example

#Intended to be run from nep-pack/ directory or nep-pack/test directory
# workspace()
# push!(LOAD_PATH, string(@__DIR__, "/../src"))
# push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
# push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))
#
# using NEPCore
# using NEPTypes
# using LinSolvers
# using NEPSolver
# using Gallery
# using IterativeSolvers
# using Base.Test

import NEPSolver.iar_chebyshev;
include("../src/method_iar_chebyshev.jl");

#explicit import needed for overloading functions from packages
import NEPCore.compute_Mlincomb



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


dep=nep_gallery("dep0");
n=size(dep,1);

IAR=@testset "IAR" begin
    @testset "accuracy eigenpairs" begin
        (λ,Q)=iar_chebyshev(dep,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100);
        @testset "IAR eigval[$i]" for i in 1:length(λ)
            @test norm(compute_Mlincomb(dep,λ[i],Q[:,i]))<eps()*100;
        end
    end

    @testset "orthogonalization" begin

    # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS
    @testset "DGKS" begin
        (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=DGKS,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100)
        @test norm(V'*V-eye(size(V,2)))<1e-6
     end

     @testset "User provided doubleGS" begin
         (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=DoubleGS,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100)
         @test norm(V'*V-eye(size(V,2)))<1e-6
      end

      @testset "ModifiedGramSchmidt" begin
          (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=ModifiedGramSchmidt,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100)
          @test norm(V'*V-eye(size(V,2)))<1e-6
      end

       @testset "ClassicalGramSchmidt" begin
           (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=ClassicalGramSchmidt,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100)
           @test norm(V'*V-eye(size(V,2)))<1e-6
       end
    end

    # Other types
    @testset "IAR CHEB PEP" begin
        srand(0); n=100; d=3;
        A = Array{Array{Float64}}(d+1)
        for j=0:d
            A[j+1]=rand(n,n)
        end
        nep=PEP(A)
        (λ,Q)=iar_chebyshev(nep,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100)

        @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
    end

    @testset "IAR CHEB GENERIC Y0" begin

        n=1000; I=[1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]; # sparsity pattern of tridiag matrix
        A0=sparse(I, J, rand(3*n-2)); A1=sparse(I, J, rand(3*n-2))
        nep=SPMF_NEP([eye(n), A0, A1],[λ->-λ,λ->eye(λ),λ->expm(-λ)])

        compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
        (λ,Q,err)=iar_chebyshev(nep,σ=0,γ=1,Neig=7,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=1,displaylevel=0,a=-1,b=2)

        @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
    end

    @testset "IAR CHEB DEP DELAYS>1" begin
        srand(0)
        n=100; A1=rand(n,n); A2=rand(n,n); A3=rand(n,n);
        tau1=0; tau2=2.3; tau3=.1;
        nep=DEP([A1,A2,A3],[tau1,tau2,tau3])
        (λ,Q)=iar_chebyshev(nep,σ=0,Neig=5,displaylevel=0,maxit=90,tol=eps()*100)

        @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
    end

    @testset "IAR CHEB DEP SHIFTED AND SCALED" begin
        nep=nep_gallery("dep0_tridiag",1000)
        (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;Neig=5,displaylevel=0,maxit=100,tol=eps()*100)

        @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
    end

    @testset "IAR CHEB PEP SHIFTED AND SCALED" begin
        srand(0)
        n=100;
        d=3;
        A = Array{Array{Float64}}(d+1)
        for j=0:d
            A[j+1]=rand(n,n)
        end
        nep=PEP(A)
        (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;Neig=5,displaylevel=0,maxit=100,tol=eps()*100)

        @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
    end

end

#TODO: add test for compute_y0 provided as input (dep+qep)
