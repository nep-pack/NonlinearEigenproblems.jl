#  Tests for the projected NEPs

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random

#using Winston # For plotting

@testset "Projected problems" begin

    @bench @testset "Problem $nepstr" for nepstr in ("pep", "dep", "sqrtm")

        local nep::NEP
        if (nepstr == "pep")
            nep=nep_gallery("pep0_sparse_003");
        elseif (nepstr == "dep")

            n=5;
            Random.seed!(1)
            A0=randn(5,5);
            A1=randn(5,5);
            t=3.0

            minusop= S-> -S
            oneop= S -> Matrix{eltype(S)}(I, size(S))
            expmop= S -> exp(Matrix(-t*S))
            fi=[minusop, oneop, expmop];

            nep = SPMF_NEP([Matrix(1.0I, n, n), A0, A1], fi)

        elseif (nepstr == "sqrtm")

            n=8;
            Random.seed!(1)
            A0=randn(n,n);
            A1=randn(n,n);
            A2=randn(n,n)/10;
            t=3.0

            minusop= S-> -S
            oneop= S -> Matrix{eltype(S)}(I, size(S))
            expmop= S -> sqrt(Matrix(-t*S) + 30*I)
            fi=[minusop, oneop, expmop];

            nep=SPMF_NEP([A0,A1,A2],fi)
        end

        #

        n=size(nep,1);
        @info "Running Newton Raphson"
        tol=1e-12
        λ,x =newton(nep,maxit=30,λ=1+1im,tol=tol,
                    v=ones(n));
        # Check residual is small
        λ_exact=λ
        @test norm(compute_Mlincomb(nep,λ,x))<tol*1000;

        # Create a projected NEP
        pnep=create_proj_NEP(nep);
        Random.seed!(1);
        V=randn(size(nep,1),2)
        Q,R=qr(hcat(V,x)) # Make the eigenspace a part of the projection subspace
        Q = Matrix(Q)
        set_projectmatrices!(pnep,Q,Q);
        @info "Running Newton on projected problem with very good start value"
        λ1,z1=newton(pnep,λ=(λ_exact+0.00001),displaylevel=0,v=ones(size(pnep,1)))

        x1=Q*z1; x1=x1/x1[1];

        # should be small since the eigenvector is in the subspace
        @test norm(x/x[1]-x1)<1e-8

        @info "Running IAR for projected problem (no special starting vector)"
        λv,X=iar(pnep,σ=complex(round(λ_exact*10)/10),displaylevel=1,
                 Neig=3,maxit=100, v=ones(size(pnep,1)))

        diff = minimum(abs.(λv .- λ_exact))
        @info "Difference of solution from projected problem: $diff"
        @test diff < sqrt(eps())
    end
end
