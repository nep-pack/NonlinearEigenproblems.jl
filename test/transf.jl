# Run tests for the NEP transformations

using NonlinearEigenproblems
using Test
using LinearAlgebra

@testset "transformations" begin

    @bench @testset "Shift and scale" begin
        orgnep=nep_gallery("dep0");

        # Test the shift by solving the orgnep and
        # check with residual of transformed nep.
        σ=-0.3+0.1im;
        nep1=shift_and_scale(orgnep,shift=σ);
        orgλ,orgv=augnewton(orgnep,v=ones(size(orgnep,1)))
        @test norm(compute_Mlincomb(nep1,orgλ-σ,orgv))<100*eps()
        λ1,v1=quasinewton(nep1,v=ones(size(orgnep,1)))
        @test abs(λ1+σ-orgλ)<eps()*100 # check that we get the same eigvals


        # Test shift and scaling
        σ=-0.4+0.01im; α=0.5
        nep2=shift_and_scale(orgnep,shift=σ,scale=α);
        λ2,v2=quasinewton(nep2,v=ones(size(nep2,1)))
        @test abs((α*λ2+σ)-orgλ)<eps()*1000


        # Check that PEP transformations correctly transform coefficients
        pep0=nep_gallery("pep0",10);
        σ=1+0.3im;
        α=3;
        pep1=shift_and_scale(pep0,shift=σ,scale=α)
        λ,v= quasinewton(pep0,λ=-0.3-0.7im,v=ones(size(pep0,1)));
        norm(compute_Mlincomb(pep0, λ, v))
        λv,V=polyeig(pep1);
        @test minimum(abs.(λv .- (λ-σ)/α)) < eps()*1000

        # Check that real PEP with real transformation is still real
        σ=3;
        α=1
        pep2=shift_and_scale(pep0,shift=σ,scale=α)
        @test !isa(Complex,eltype(pep2.A[1])) # Preserve realness


        nep3=nep_gallery("qdep0")
        λ,v= quasinewton(nep3,λ=1+1im,v=ones(size(nep3,1)));
        σ=-3+0.3im
        α=0.9;
        nep3_transf=shift_and_scale(nep3,shift=σ,scale=α);
        @test norm(compute_Mlincomb(nep3_transf,(λ-σ)/α,v))<sqrt(eps())*1000;
        λ,V=iar(nep3_transf,σ=0,neigs=2,maxit=60,v=ones(size(nep3_transf,1)))
        @test norm(compute_Mlincomb(nep3,α*λ[1]+σ,V[:,1]))<sqrt(eps())
        @test norm(compute_Mlincomb(nep3,α*λ[2]+σ,V[:,2]))<sqrt(eps())


        λ,V=iar(nep3_transf,σ=0,neigs=2,maxit=60,v=ones(size(nep3_transf,1)))
        @test norm(compute_Mlincomb(nep3,α*λ[1]+σ,V[:,1]))<sqrt(eps())
        @test norm(compute_Mlincomb(nep3,α*λ[2]+σ,V[:,2]))<sqrt(eps())
        #

    end
    @bench @testset "Möbius transformation" begin
        # Checks mobius_transformation

        pep0=nep_gallery("pep0")
        a=3+0.1im
        b=0.1im;
        c=1;
        d=1-0.3im
        pep0_transf=mobius_transform(pep0,a=a,b=b,c=c,d=d)
        λ,v= quasinewton(pep0_transf,λ=0.3im,v=ones(size(pep0,1)));
        λorg=(a*λ+b)/(c*λ+d)
        @test norm(compute_Mlincomb(pep0,λorg,v))<sqrt(eps());


        nep4=nep_gallery("qdep0")
        a=3-0.01im
        b=0.1im;
        c=1;
        d=1-0.3im
        nep4_transf=mobius_transform(nep4,a=a,b=b,c=c,d=d)
        @test isa(nep4_transf,SPMF_NEP)  # Möbius transformed SPMF-NEP is still a SPMF
        λ,v= quasinewton(nep4_transf,λ=1-0.3im,v=ones(size(nep4,1)));
        λorg=(a*λ+b)/(c*λ+d)
        @test norm(compute_Mlincomb(nep4,λorg,v))<sqrt(eps())*10;


        n=size(nep4,1);
        V=(1:n)*(1:5)'/n; # Arbitrary matrices
        S=(1:5)*(1:5)'/5-3I;
        s = Matrix(3.0*I, 1, 1)
        W1 = compute_MM(nep4, (c*S + d*I) \ (a*S + b*I), V)
        W2=compute_MM(nep4_transf,S,V);
        @test opnorm(W1-W2)<sqrt(eps())

        # Check that mobius_transformation becomes shift_and_scale
        c=0; d=1;
        nep4_transf1=mobius_transform(nep4,a=a,b=b,c=c,d=d)
        nep4_transf2=shift_and_scale(nep4,scale=a,shift=b);
        λ=3;
        w=compute_Mlincomb(nep4,λ*a+b,V[:,1])
        w1=compute_Mlincomb(nep4_transf1,λ,V[:,1])
        w2=compute_Mlincomb(nep4_transf2,λ,V[:,1])
        @test norm(w-w1)<eps()*100
        @test norm(w-w2)<eps()*100
        @test norm(w1-w2)<eps()*100
        W1=compute_Mlincomb(nep4_transf1,λ,V,[0,1,0,1,0])
        W2=compute_Mlincomb(nep4_transf2,λ,V,[0,1,0,1,0]);
        @test norm(W1-W2)<sqrt(eps())





    end

    @bench @testset "taylor exp" begin
        nep=nep_gallery("dep0");
        d=5;
        pep=taylor_expansion_pep(nep,d);
        s=0.1;
        @test norm(compute_Mder(pep,s)-compute_Mder(nep,s)) < abs(s)^d
    end



end
