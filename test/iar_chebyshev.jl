using NonlinearEigenproblems
using Test
using Random
using SparseArrays
using LinearAlgebra
using IterativeSolvers

#explicit import needed for overloading functions from packages
import NonlinearEigenproblems.NEPCore.compute_Mlincomb

Random.seed!(0);

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
function orthogonalize_and_normalize!(V,v,h,::DoubleGS)
    doubleGS_function!(V, v, h) end


# The user can provide compute_y0_cheb and associated procumputation
# Here there is a naive example for the QEP. Relate to the test "compute_y0 AS INPUT FOR QEP (naive)"
# Import and create a type to overload the compute_y0_cheb function
import NonlinearEigenproblems.NEPSolver.ComputeY0Cheb
import NonlinearEigenproblems.NEPSolver.AbstractPrecomputeData
import NonlinearEigenproblems.NEPSolver.precompute_data
abstract type ComputeY0Cheb_QEP <: NEPSolver.ComputeY0Cheb end
struct PrecomputeData_QEP <: AbstractPrecomputeData
    precomp_PEP; nep_pep
end


function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QEP},a,b,m,γ,σ)
    A0,A1,A2=get_Av(nep)
    nep_pep=PEP([A0,A1,A2])
    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);
    return PrecomputeData_QEP(precomp_PEP,nep_pep)
end
# overload the function compute_y0_cheb
import NonlinearEigenproblems.NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QEP},x,y,M0inv,precomp::PrecomputeData_QEP)
    return compute_y0_cheb(T,precomp.nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)
end

# A less trivial example on how the user can provide his own compute_y0_cheb and associated procumputation
# Here there is a naive example for the QDEP. Relate to the test "compute_y0 AS INPUT FOR QEP (naive)"
abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
struct PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_PEP; precomp_DEP; nep_pep; nep_dep
end
function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (I + λ I - λ^2 I) + (-λ A0 + A1 exp(-τ λ))
    A0,A1=get_Av(nep)[2:3]
    n=size(nep,1)
    J = Matrix{T}(I,n,n)
    nep_pep=PEP([J, J, -J]) # the PEP part is defined as
    nep_dep=DEP([A0,A1],[0.0,1.0]);     # the DEP part is defined as
    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);
    precomp_DEP=precompute_data(T,nep_dep,NEPSolver.ComputeY0ChebDEP,a,b,m,γ,σ);
    # combine the precomputations
    return PrecomputeData_QDEP(precomp_PEP,precomp_DEP,nep_pep,nep_dep)
end
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    return compute_y0_cheb(T,precomp.nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)+ compute_y0_cheb(T,precomp.nep_dep,NEPSolver.ComputeY0ChebDEP,x,y,M0inv,precomp.precomp_DEP)+y*(view(precomp.precomp_DEP.Tc,1:size(x,2)+1));
    # y*Tc subtracted at each call of compute_y0iar_cheb. Therefore since we do two calls, we need to add it back once.
end

@testset "IAR Chebyshev version" begin
    dep=nep_gallery("neuron0"); n=size(dep,1);
    @bench @testset "Scale Cheb's to different interval w DEP" begin
        a=-maximum(dep.tauv)
        b=0;
        (λ,Q)=iar_chebyshev(dep,a=a,b=b,neigs=10,maxit=100,v=ones(size(dep,1)))
        verify_lambdas(10, dep, λ, Q, n*sqrt(eps()))
    end


    dep=nep_gallery("dep0"); n=size(dep,1);
    @bench @testset "accuracy eigenpairs" begin
        (λ,Q)=iar_chebyshev(dep,σ=0,neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(dep,1)));
        verify_lambdas(5, dep, λ, Q, n*sqrt(eps()))
    end

    @testset "Compute as many eigenpairs as possible (neigs=Inf)" begin
        (λ,Q)=iar_chebyshev(dep,σ=0,neigs=Inf,logger=0,maxit=30,
                            tol=eps()*100,v=ones(size(dep,1)));
        verify_lambdas(8, dep, λ, Q, n*sqrt(eps()))
    end



    @testset "orthogonalization" begin

    # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS
    @bench @testset "DGKS" begin
        (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=DGKS(),σ=0,neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(dep,1)))
        @test opnorm(V'*V - Matrix(1.0I, size(V,2), size(V,2))) < n*sqrt(eps())
     end

     @bench @testset "User provided doubleGS" begin
         (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=DoubleGS(),σ=0,neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(dep,1)))
         @test opnorm(V'*V - Matrix(1.0I, size(V,2), size(V,2))) < n*sqrt(eps())
      end

      @bench @testset "ModifiedGramSchmidt" begin
          (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=ModifiedGramSchmidt(),σ=0,neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(dep,1)))
          @test opnorm(V'*V - Matrix(1.0I, size(V,2), size(V,2))) < n*sqrt(eps())
      end

       @bench @testset "ClassicalGramSchmidt" begin
           (λ,Q,err,V)=iar_chebyshev(dep,orthmethod=ClassicalGramSchmidt(),σ=0,neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(dep,1)))
           @test opnorm(V'*V - Matrix(1.0I, size(V,2), size(V,2))) < n*sqrt(eps())
       end
    end

    # Other types
    @testset "compute_y0_method for different types" begin
        @bench @testset "PEP" begin
            n=100; d=3;
            A = Array{Array{Float64}}(undef, d+1)
            for j=0:d
                A[j+1]=rand(n,n)
            end
            nep=PEP(A)
            (λ,Q)=iar_chebyshev(nep,σ=0,neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(nep,1)))

            verify_lambdas(5, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "DEP WITH GENERIC Y0" begin
            n=1000; K=[1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]; # sparsity pattern of tridiag matrix
            A0=sparse(K, J, rand(3*n-2)); A1=sparse(K, J, rand(3*n-2))
            nep = SPMF_NEP([sparse(1.0I, n, n), A0, A1], [λ -> -λ, λ -> one(λ), λ -> exp(-λ)])

            compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
            (λ,Q,err)=iar_chebyshev(nep,σ=0,γ=1,neigs=7,logger=0,maxit=100,tol=eps()*100,check_error_every=1,a=-1,b=2,v=ones(size(nep,1)))

            verify_lambdas(7, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "DEP WITH DELAYS>1" begin
            n=100; A1=rand(n,n); A2=rand(n,n); A3=rand(n,n);
            tau1=0; tau2=1.1; tau3=.2;
            nep=DEP([A1,A2,A3],[tau1,tau2,tau3])
            (λ,Q)=iar_chebyshev(nep,σ=0,neigs=3,logger=0,maxit=90,tol=eps()*100,v=ones(size(nep,1)))

            verify_lambdas(3, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "DEP SHIFTED AND SCALED" begin
            nep=nep_gallery("dep0_tridiag",1000)
            (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(nep,1)))
            verify_lambdas(5, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "PEP SHIFTED AND SCALED" begin
            n=100; d=3;
            A = Array{Array{Float64}}(undef, d+1)
            for j=0:d
                A[j+1]=rand(n,n)
            end
            nep=PEP(A)
            (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;neigs=5,logger=0,maxit=100,tol=eps()*100,v=ones(size(nep,1)))

            verify_lambdas(5, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "compute_y0 AS INPUT FOR QEP (naive)" begin
            A0=rand(n,n); A1=rand(n,n); A2=rand(n,n);
            nep = SPMF_NEP([A0, A1, A2], [λ -> one(λ), λ -> λ, λ -> λ^2])
            λ,Q,err,V = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QEP,maxit=100,neigs=10,σ=0.0,γ=1,logger=0,check_error_every=1,v=ones(size(nep,1)));

            verify_lambdas(10, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "QDEP IN SPMF format" begin
            nep=nep_gallery("qdep1")
            λ,Q,err,V = iar_chebyshev(nep,maxit=100,neigs=8,σ=0.0,γ=1,logger=0,check_error_every=1,v=ones(size(nep,1)));
            verify_lambdas(8, nep, λ, Q, n*sqrt(eps()))
            # with scaling
            λ,Q,err,V = iar_chebyshev(nep,maxit=100,neigs=8,σ=0.0,γ=0.9,logger=0,check_error_every=1,v=ones(size(nep,1)));
            verify_lambdas(8, nep, λ, Q, n*sqrt(eps()))
        end

        @bench @testset "PEP in SPMF format" begin
            A0=rand(n,n); A1=rand(n,n); A2=rand(n,n);
            nep = SPMF_NEP([A0, A1, A2], [λ -> one(λ), λ -> λ, λ -> λ^2])

            λ,Q,err,V = iar_chebyshev(nep,maxit=100,neigs=10,σ=0.0,γ=1,logger=0,check_error_every=1,v=ones(n));
            verify_lambdas(10, nep, λ, Q, n*sqrt(eps()))

            λ2,Q2,err2,V2, H2 = iar_chebyshev(nep,maxit=100,neigs=20,σ=0.0,γ=1,logger=0,check_error_every=1,compute_y0_method=ComputeY0Cheb,v=ones(n));
            @test opnorm(V[:,1:10]-V2[:,1:10])<n*sqrt(eps());
        end

        @bench @testset "DEP format with ComputeY0ChebSPMF_NEP" begin
            nep=nep_gallery("dep0_tridiag",1000)
            (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2;neigs=5,logger=0,maxit=100,tol=eps()*100,compute_y0_method=NEPSolver.ComputeY0ChebSPMF_NEP)
            verify_lambdas(5, nep, λ, Q, 1e-10)
        end

        @bench @testset "PEP format with ComputeY0ChebSPMF_NEP" begin
            n=100; d=3;
            A = Array{Array{Float64}}(undef, d+1)
            Random.seed!(0)
            for j=0:d
                A[j+1]=rand(n,n)
            end
            nep=PEP(A)
            (λ,Q)=iar_chebyshev(nep,σ=-1,γ=2,neigs=5,logger=0,maxit=100,tol=eps()*100,compute_y0_method=NEPSolver.ComputeY0ChebSPMF_NEP,v=ones(size(nep,1)))

            verify_lambdas(5, nep, λ, Q, 1e-10)
        end

        @bench @testset "compute_y0 AS INPUT FOR QDEP (combine DEP and PEP)" begin
            nep=nep_gallery("qdep1")
            λ,Q,err,V,H = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=100,neigs=10,logger=0,
                                        check_error_every=1,v=ones(size(nep,1)));
            verify_lambdas(10, nep, λ, Q, 1e-10)

            λ2,Q2,err2,V2,H2 = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=100,neigs=10,logger=0,
                                        check_error_every=1,v=ones(size(nep,1)));
            verify_lambdas(10, nep, λ, Q, 1e-10)

            @test opnorm(V[:,1:10]-V2[:,1:10])+opnorm(H[1:10,1:10]-H2[1:10,1:10])<1e-10;
        end
    end

    @testset "Errors thrown" begin
        np=100;
        dep=nep_gallery("dep0",np);
        @test_throws NEPCore.NoConvergenceException (λ,Q)=iar_chebyshev(dep,σ=0,neigs=8,logger=0,maxit=10,tol=eps()*100,v=ones(size(dep,1)));
    end
end
