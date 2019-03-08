# Run tests for the deflation

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra

@bench @testset "Deflation (combined with MSLP)" begin

    nep=nep_gallery("dep0");
    nep = DEP(nep.n, nep.A, [0, 0.8])

    n=size(nep,1);

    (λ,v)=newton(nep,v=ones(n))

    v=v/norm(v);
    dnep=deflate_eigpair(nep,λ,v)

    for k=1:3
        λ0=-0.1+0.1im;

        (λ2,v2)=mslp(dnep,maxit=1000,
                     λ=λ0,
                     displaylevel=0
                     );
        @info "Computing eigenvalue k=$k"
        v2=v2/norm(v2);
        # Expand the partial Schur factorization with the computed solution
        dnep=deflate_eigpair(dnep,λ2,v2);
    end

    (λv,V)=get_deflated_eigpairs(dnep);

    # Test that λv and V are now eigpairs
    for i=1:size(λv,1)
        v=V[:,i]/norm(V[:,i]);
        λ=λv[i];
        @test norm(compute_Mlincomb(nep,λ,v))<sqrt(eps())
    end
end


@bench @testset "Deflation modes (new)" begin

    nep=nep_gallery("dep0_sparse",3);
    n=size(nep,1);
    local λ,v;
    (λ,v)=augnewton(nep,v=ones(n),λ=1.65,tol=1e-11)
    dnep=deflate_eigpair(nep,λ,v,mode=:Generic);
    dnep2=deflate_eigpair(nep,λ,v,mode=:SPMF);
    dnep3=deflate_eigpair(nep,λ,v,mode=:MM);
    λ=2+2im;
    @test compute_Mder(dnep,λ) ≈ compute_Mder(dnep2,λ)
    @test compute_Mder(dnep,λ) ≈ compute_Mder(dnep3,λ)

    X=[1 2; 3 4 ; 5 5.0; -1 -1];

    @test compute_Mlincomb(dnep,λ,X) ≈ compute_Mlincomb(dnep3,λ,X)
    @test compute_Mlincomb(dnep,λ,X) ≈ compute_Mlincomb(dnep2,λ,X)

end
