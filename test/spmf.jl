# Tests for SPMF-code

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random
using SparseArrays

@testset "SPMF" begin
    # Create an SPMF
    n=5;
    Random.seed!(1)
    A0=sparse(randn(n,n));
    A1=sparse(randn(n,n));
    t=3.0

    minusop = S -> -S
    oneop = S -> one(S)
    expmop = S -> exp(-t*S)
    fi=[minusop, oneop, expmop];

    J = SparseMatrixCSC(1.0I, n, n)
    nep1 = SPMF_NEP([J, A0, A1], fi)
    nep2 = SPMF_NEP([J, A0, A1], fi; Schur_fact=true)

    @bench @testset "compute_MM" begin

        S=randn(3,3);
        V=randn(n,3);

        # Check if prefactorize with Schur gives the same result
        @test opnorm(compute_MM(nep1,S,V)-compute_MM(nep2,S,V))<sqrt(eps())
        # Check compute_MM
        @test opnorm(compute_MM(nep1,S,V)-(-V*S+A0*V+A1*V*exp(-t*S)))<sqrt(eps())


        # Check if compute_MM is correct (by comparing against diagonalization of S).
        # This is done with the identity MM(V,S)=MM(V,W*D*inv(W))=MM(V*W,D)*inv(W)
        # and MM(V1,D) for diagonal D can be computed with compute_Mlincomb

        N1=compute_MM(nep1,S,V);
        # Diagonalize S
        d,W = eigen(S)
        D = diagm(0 => d)
        V1 = V*W
        #
        N2=hcat(compute_Mlincomb(nep1,d[1],V1[:,1]),
                compute_Mlincomb(nep1,d[2],V1[:,2]),
                compute_Mlincomb(nep1,d[3],V1[:,3]))*inv(W)
        @test opnorm(N1-N2)<sqrt(eps())

    end
    @bench @testset "compute_Mder_from_MM" begin

        S=randn(3,3);
        V=randn(n,3);

        # Check compute_Mder_from_MM()
        λ=2
        T1=compute_Mder_from_MM(nep1,λ,1)
        T2 = -I - t*A1*exp(-t*λ)
        @test opnorm(T1-T2)<sqrt(eps())

        λ=2
        T3=compute_Mder_from_MM(nep1,λ,2)
        T4=t^2*A1*exp(-t*λ)
        @test opnorm(T3-T4)<sqrt(eps())



        # Check consistency of MM and Mder

        # Exact Mlincomb:
        Zexact = (-λ*I + A0 + A1*exp(-t*λ)) * V[:,1] +
            (-I - t*A1*exp(-t*λ)) * V[:,2] +
            (t^2*A1*exp(-t*λ)) * V[:,3]

        Z1=compute_Mlincomb_from_MM(nep1,λ,V,[1.0,1.0,1.0])
        @test norm(Z1-Zexact)<sqrt(eps())

        Z2=compute_Mlincomb_from_Mder(nep1,λ,V,[1.0,1.0,1.0])
        @test norm(Z2-Zexact)<sqrt(eps())

    end

    @bench @testset "Compute_Mlincomb" begin

        Random.seed!(10)
        S=randn(3,3)+120^2*I;
        V=randn(n,3);

        # Same nonlinearities as the GUN NLEVP problem
        minusop= S-> -S
        oneop = S -> Matrix(1.0I, size(S))
        sqrt1op = S -> 1im * sqrt(Matrix(S))
        sqrt2op = S -> 1im * sqrt(Matrix(S) - 108.8774^2*I)

        A0=sparse(randn(n,n))+1im*sparse(randn(n,n));
        A1=sparse(randn(n,n))+1im*sparse(randn(n,n));
        A2=sparse(randn(n,n))+1im*sparse(randn(n,n));
        A3=sparse(randn(n,n))+1im*sparse(randn(n,n));
        nep2=SPMF_NEP([A0,A1,A2,A3],[oneop,minusop,sqrt1op,sqrt2op]);

        N1=compute_MM(nep2,S,V);
        # Diagonalize S
        d,W = eigen(S)
        D = diagm(0 => d)
        V1 = V*W
        #
        N2=hcat(compute_Mlincomb(nep2,d[1],V1[:,1]),
                compute_Mlincomb(nep2,d[2],V1[:,2]),
                compute_Mlincomb(nep2,d[3],V1[:,3]))*inv(W)
        @test norm(N1-N2)<sqrt(eps())*100

    end

    @bench @testset "PEP" begin
        Random.seed!(99)
        A0=randn(5,5)
        A1=randn(5,5)
        A2=randn(5,5)
        A3=randn(5,5)
        # Check that PEP generates the
        # same values as an SPMF generated from the PEP
        nep0=PEP([A0,A1,A2, A3])
        nep1=SPMF_NEP(get_Av(nep0),get_fv(nep0))
        for  λ in [3,10,-3,-100]
            M0=compute_Mder(nep0,λ)
            M1=compute_Mder(nep1,λ)
            @test M0 ≈ M1

            M0=compute_Mder(nep0,λ,2)
            M1 = 2 .* A2 + (λ*6) .* A3
            @test M0 ≈ M1
        end

    end

    @bench @testset "REP" begin
        Random.seed!(10)
        A0=randn(5,5)
        A1=randn(5,5)
        A2=randn(5,5)
        # Check that REP generates the
        # same values as an SPMF generated from the REP
        nep0=REP([A0,A1,A2],[1,2,3.3])
        nep1=SPMF_NEP(get_Av(nep0),get_fv(nep0))
        for  λ in [3,10,-3,-100]
            M0=compute_Mder(nep0,λ)
            M1=compute_Mder(nep1,λ)
            @test M0 ≈ M1
        end

    end

    @bench @testset "DEP" begin
        Random.seed!(88)
        A1=randn(5,5)
        A2=randn(5,5)
        τ1 = 1.5
        τ2 = 3.75
        # Check that DEP generates the
        # same values as an SPMF generated from the DEP
        nep0=DEP([A1,A2],[τ1, τ2])
        nep1=SPMF_NEP(get_Av(nep0),get_fv(nep0))
        for  λ in [3,10,-3,-100]
            M0=compute_Mder(nep0,λ)
            M1=compute_Mder(nep1,λ)
            @test M0 ≈ M1
            M2 = -λ*I + A1 .* exp(-τ1*λ) + A2 .* exp(-τ2*λ)
            @test M0 ≈ M2
        end

    end

end
