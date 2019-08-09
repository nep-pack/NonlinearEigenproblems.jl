#  Tests for the projected NEPs

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random

#using Winston # For plotting

@testset "Projected problems" begin

    @bench @testset "Problem $nepstr" for nepstr in ("pep", "dep", "sqrtm")

        local nep::NEP
        if (nepstr == "pep")
            nep=nep_gallery("pep0_sparse");
        elseif (nepstr == "dep")

            n=5;
            Random.seed!(1)
            A0=randn(5,5);
            A1=randn(5,5);
            t::Float64=3.0

            minusop= S-> -S
            oneop= S -> one(S)
            expmop= S -> exp(-t*S)
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
            oneop= S -> one(S)  # one(S) replaced by I it would not return a matrix
            expmop= S -> sqrt(-t*S + 30*one(S))
            fi=[minusop, oneop, expmop];

            nep=SPMF_NEP([A0,A1,A2],fi)
        end

        #

        n=size(nep,1);
        @info "Running Newton Raphson ($nepstr)"
        tol=1e-12
        λ,x =newton(nep,maxit=30,λ=1+1im,tol=tol,
                    v=ones(n));
        # Check residual is small
        λ_exact=λ
        @test norm(compute_Mlincomb(nep,λ,x))<tol*1000;

        # Create a projected NEP
        pnep=create_proj_NEP(nep,4); # maxsize=4
        Random.seed!(1);
        V=randn(size(nep,1),2)
        Q,R=qr(hcat(V,x)) # Make the eigenspace a part of the projection subspace
        Q = Matrix(Q)
        set_projectmatrices!(pnep,Q,Q);
        @info "Running Newton on projected problem with very good start value ($nepstr)"
        λ1,z1=newton(pnep,λ=(λ_exact+0.00001),logger=0,v=ones(size(pnep,1)))

        x1=Q*z1; x1=x1/x1[1];



        # should be small since the eigenvector is in the subspace
        @test norm(x/x[1]-x1)<1e-8

        @info "Running IAR for projected problem ($nepstr)"
        λv,X=iar(pnep,σ=complex(round(λ_exact*10)/10),
                 neigs=3,maxit=100, v=ones(size(pnep,1)))

        mydiff = minimum(abs.(λv .- λ_exact))
        @info "Difference of solution from projected problem ($nepstr): $mydiff"
        @test mydiff < sqrt(eps())


        @info "Testing to expand and running Newton ($nepstr)"
        n=size(nep,1);
        Vnew=[Q ones(n)];
        Wnew=[Q ones(n)];
        expand_projectmatrices!(pnep,Wnew,Vnew);
        λ1,zz1=newton(pnep,λ=(λ_exact+0.00001),logger=0,
                      v=ones(4),maxit=20)

        xx1=Vnew*zz1; xx1=xx1/xx1[1]

        # should be small since the eigenvector is in the subspace
        @test norm(x/x[1]-x1)<1e-8

    end
end
