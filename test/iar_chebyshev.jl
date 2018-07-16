#Run tests for the dep_distributed example

#Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using IterativeSolvers
using Base.Test

#import NEPSolver.iar_chebyshev;
#include("../src/method_iar_chebyshev.jl");

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


# The user can provide compute_y0_cheb and associated procumputation
# Here there is a naive example for the QEP. Relate to the test "compute_y0 AS INPUT FOR QEP (naive)"
# Import and create a type to overload the compute_y0_cheb function
import NEPSolver.ComputeY0Cheb
import NEPSolver.AbstractPrecomputeData
import NEPSolver.precompute_data
abstract type ComputeY0Cheb_QEP <: NEPSolver.ComputeY0Cheb end
type PrecomputeData_QEP <: AbstractPrecomputeData
    precomp_PEP; nep_pep
end
function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QEP},a,b,m,γ,σ)
    A0,A1,A2=get_Av(nep)
    nep_pep=PEP([A0,A1,A2])
    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);
    return PrecomputeData_QEP(precomp_PEP,nep_pep)
end
# overload the function compute_y0_cheb
import NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QEP},x,y,M0inv,precomp::PrecomputeData_QEP)
    return compute_y0_cheb(T,precomp.nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)
end

# A less trivial example on how the user can provide his own compute_y0_cheb and associated procumputation
# Here there is a naive example for the QDEP. Relate to the test "compute_y0 AS INPUT FOR QEP (naive)"
abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
type PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_PEP; precomp_DEP; nep_pep; nep_dep
end
function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (I + λ I - λ^2 I) + (-λ A0 + A1 exp(-τ λ))
    A0,A1=get_Av(nep)[2:3]
    n=size(nep,1)
    nep_pep=PEP([eye(T,n,n), eye(T,n,n), -eye(T,n,n)]); # the PEP part is defined as
    nep_dep=DEP([A0,A1],[0.0,1.0]);     # the DEP part is defined as

    # NOTE: it should be enough this in order to pass the test
    if σ!=zero(T) || γ!=one(T)
        nep_pep=shift_and_scale(nep_pep,shift=σ,scale=γ)
        nep_dep=shift_and_scale(nep_dep,shift=σ,scale=γ)
        σ=zero(T); γ=one(T)
    end

    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);
    precomp_DEP=precompute_data(T,nep_dep,NEPSolver.ComputeY0ChebDEP,a,b,m,γ,σ);
    # combine the precomputations
    return PrecomputeData_QDEP(precomp_PEP,precomp_DEP,nep_pep,nep_dep)
end
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    return compute_y0_cheb(T,precomp.nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)+ compute_y0_cheb(T,precomp.nep_dep,NEPSolver.ComputeY0ChebDEP,x,y,M0inv,precomp.precomp_DEP)+y*(view(precomp.precomp_DEP.Tc,1:size(x,2)+1));
    # y*Tc subtracted at each call of compute_y0iar_cheb. Therefore since we do two calls, we need to add it back once.
end


IAR_cheb=@testset "IAR Chebyshev version" begin
    dep=nep_gallery("dep0"); n=size(dep,1);
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
    @testset "compute_y0_method for different types" begin

        @testset "PEP" begin
            srand(0); n=100; d=3;
            A = Array{Array{Float64}}(d+1)
            for j=0:d
                A[j+1]=rand(n,n)
            end
            nep=PEP(A)
            (λ,Q)=iar_chebyshev(nep,σ=0,Neig=5,displaylevel=0,maxit=100,tol=eps()*100)

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end

        @testset "DEP WITH GENERIC Y0" begin

            n=1000; I=[1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]; # sparsity pattern of tridiag matrix
            A0=sparse(I, J, rand(3*n-2)); A1=sparse(I, J, rand(3*n-2))
            nep=SPMF_NEP([eye(n), A0, A1],[λ->-λ,λ->eye(λ),λ->expm(-λ)])

            compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
            (λ,Q,err)=iar_chebyshev(nep,σ=0,γ=1,Neig=7,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=1,displaylevel=0,a=-1,b=2)

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end

        @testset "DEP WITH DELAYS>1" begin
            srand(0)
            n=100; A1=rand(n,n); A2=rand(n,n); A3=rand(n,n);
            tau1=0; tau2=2.3; tau3=.1;
            nep=DEP([A1,A2,A3],[tau1,tau2,tau3])
            (λ,Q)=iar_chebyshev(nep,σ=0,Neig=5,displaylevel=0,maxit=90,tol=eps()*100)

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end

        @testset "DEP SHIFTED AND SCALED" begin
            nep=nep_gallery("dep0_tridiag",1000)
            (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;Neig=5,displaylevel=0,maxit=100,tol=eps()*100)

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end

        @testset "PEP SHIFTED AND SCALED" begin
            srand(0); n=100; d=3;
            A = Array{Array{Float64}}(d+1)
            for j=0:d
                A[j+1]=rand(n,n)
            end
            nep=PEP(A);       (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;Neig=5,displaylevel=0,maxit=100,tol=eps()*100)

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end

        @testset "compute_y0 AS INPUT FOR QEP (naive)" begin
            srand(0);   A0=rand(n,n); A1=rand(n,n); A2=rand(n,n);
            nep=SPMF_NEP([A0, A1, A2],[λ->eye(λ),λ->λ,λ->λ^2])

            λ,Q,err,V = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QEP,maxit=100,Neig=10,σ=0.0,γ=1,displaylevel=0,check_error_every=1);
            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end

        @testset "QDEP IN SPMF format" begin
            nep=nep_gallery("qdep1")
            λ,Q,err,V = iar_chebyshev(nep,maxit=100,Neig=8,σ=0.0,γ=1,displaylevel=0,check_error_every=1);
            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
            # with scaling
            λ,Q,err,V = iar_chebyshev(nep,maxit=100,Neig=8,σ=0.0,γ=0.9,displaylevel=0,check_error_every=1);
            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;

        end

        @testset "PEP in SPMF format" begin
            srand(0);   A0=rand(n,n); A1=rand(n,n); A2=rand(n,n);
            nep=SPMF_NEP([A0, A1, A2],[λ->eye(λ),λ->λ,λ->λ^2])

            λ,Q,err,V = iar_chebyshev(nep,maxit=100,Neig=10,σ=0.0,γ=1,displaylevel=0,check_error_every=1,v=ones(n));
            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;

            λ2,Q2,err2,V2, H2 = iar_chebyshev(nep,maxit=100,Neig=20,σ=0.0,γ=1,displaylevel=0,check_error_every=1,compute_y0_method=ComputeY0Cheb,v=ones(n));

            @test norm(V[:,1:10]-V2[:,1:10])<1e-6;

        end

        @testset "DEP format with ComputeY0ChebSPMF_NEP" begin
            nep=nep_gallery("dep0_tridiag",1000)
            (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;Neig=5,displaylevel=0,maxit=100,tol=eps()*100,compute_y0_method=NEPSolver.ComputeY0ChebSPMF_NEP)

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10
        end


        @testset "compute_y0 AS INPUT FOR QDEP (combine DEP and PEP)" begin
            nep=nep_gallery("qdep1")
            λ,Q,err,V,H = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,
                                        maxit=100,Neig=10,σ=0.0,γ=0.5,displaylevel=0,
                                        check_error_every=1,v=ones(size(nep,1)));

            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;


            λ2,Q2,err2,V2,H2 = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,
                                             maxit=100,Neig=10,σ=0.0,γ=0.5,displaylevel=0,
                                             check_error_every=1,v=ones(size(nep,1)));

            @test norm(V[:,1:10]-V2[:,1:10])<1e-6;


            # Check scaling for QDEP
            # NOTE: this test is failing because the scaling parameter γ has to be one. The user-provided method compute_y0_method does not depend on γ and it is implicitly assumed γ=1. Fix this function or change back γ=1 in order to pass this test.
            λ,Q,err,V3,H = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,
                                        maxit=100,Neig=10,σ=0.0,γ=0.9,displaylevel=0,
                                        check_error_every=1,v=ones(size(nep,1)));
            @test compute_resnorm(nep,λ[1],Q[:,1])<1e-10;
        end


    end
end



1;
