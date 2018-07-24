
# Tests for SPMF-code


# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_running_all_tests) || global_running_all_tests != true
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery

    using Base.Test
end

@testset "SPMF" begin
    # Create an SPMF
    n=5;
    srand(1)
    A0=sparse(randn(n,n));
    A1=sparse(randn(n,n));
    t=3.0

    minusop= S-> -S
    oneop= S -> eye(S)
    expmop= S -> expm(full(-t*S))
    fi=[minusop, oneop, expmop];

    nep1=SPMF_NEP([speye(n),A0,A1],fi)
    nep2=SPMF_NEP([speye(n),A0,A1],fi, true)


    @testset "compute_MM" begin

        S=randn(3,3);
        V=randn(n,3);

        # Check if prefactorize with Schur gives the same result
        @test norm(compute_MM(nep1,S,V)-compute_MM(nep2,S,V))<sqrt(eps())
        # Check compute_MM
        @test norm(compute_MM(nep1,S,V)-(-V*S+A0*V+A1*V*expm(-t*S)))<sqrt(eps())


        # Check if compute_MM is correct (by comparing against diagonalization of S).
        # This is done with the identity MM(V,S)=MM(V,W*D*inv(W))=MM(V*W,D)*inv(W)
        # and MM(V1,D) for diagonal D can be computed with compute_Mlincomb

        N1=compute_MM(nep1,S,V);
        # Diagonalize S
        (d,W)=eig(S);
        D=diagm(d);
        V1=V*W;
        #
        N2=hcat(compute_Mlincomb(nep1,d[1],V1[:,1]),
                compute_Mlincomb(nep1,d[2],V1[:,2]),
                compute_Mlincomb(nep1,d[3],V1[:,3]))*inv(W)
        @test norm(N1-N2)<sqrt(eps())

    end
    @testset "compute_Mder_from_MM" begin

        S=randn(3,3);
        V=randn(n,3);

        # Check compute_Mder_from_MM()
        λ=2
        T1=compute_Mder_from_MM(nep1,λ,1)
        T2=-1*speye(n)-t*A1*expm(-t*λ)
        @test norm(T1-T2)<sqrt(eps())

        λ=2
        T3=compute_Mder_from_MM(nep1,λ,2)
        T4=t^2*A1*expm(-t*λ)
        @test norm(T3-T4)<sqrt(eps())



        # Check consistency of MM and Mder

        # Exact Mlincomb:
        Zexact=(-λ*speye(n)+A0+A1*exp(-t*λ))*V[:,1]+
        (-speye(n)-t*A1*exp(-t*λ))*V[:,2]+
        (t^2*A1*exp(-t*λ))*V[:,3];

        Z1=compute_Mlincomb_from_MM(nep1,λ,V,[1.0,1.0,1.0])
        @test norm(Z1-Zexact)<sqrt(eps())

        Z2=compute_Mlincomb_from_Mder(nep1,λ,V,[1.0,1.0,1.0])
        @test norm(Z2-Zexact)<sqrt(eps())

    end

    @testset "Compute_Mlincomb" begin

        S=randn(3,3);
        V=randn(n,3);

        # Same nonlinearities as the GUN NLEVP problem
        minusop= S-> -S
        oneop= S -> eye(size(S,1),size(S,2))
        sqrt1op= S -> 1im*sqrtm(full(S))
        sqrt2op= S -> 1im*sqrtm(full(S)-108.8774^2*eye(S))

        A0=sparse(randn(n,n))+1im*sparse(randn(n,n));
        A1=sparse(randn(n,n))+1im*sparse(randn(n,n));
        A2=sparse(randn(n,n))+1im*sparse(randn(n,n));
        A3=sparse(randn(n,n))+1im*sparse(randn(n,n));
        nep2=SPMF_NEP([A0,A1,A2,A3],[oneop,minusop,sqrt1op,sqrt2op]);

        N1=compute_MM(nep2,S,V);
        # Diagonalize S
        (d,W)=eig(S);
        D=diagm(d);
        V1=V*W;
        #
        N2=hcat(compute_Mlincomb(nep2,d[1],V1[:,1]),
                compute_Mlincomb(nep2,d[2],V1[:,2]),
                compute_Mlincomb(nep2,d[3],V1[:,3]))*inv(W)
        @test norm(N1-N2)<sqrt(eps())*100

    end

    @testset "REP" begin
        srand(10)
        A0=randn(5,5);
        A1=randn(5,5);
        A2=randn(5,5);
        # Check that REP generates the
        # same values as an SPMF generated from the REP
        nep0=REP([A0,A1,A2],[1,2,3.3])
        nep1=SPMF_NEP(get_Av(nep0),get_fv(nep0))
        for  λ in [3,10,-3,-100];
            M0=compute_Mder(nep0,λ)
            M1=compute_Mder(nep1,λ)
            @test M0 ≈ M1
        end

    end

end
