# Run tests for the deflation

using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.NEPTypes
using NonlinearEigenproblems.Gallery
using Test
using LinearAlgebra

@testset "Deflation (combined with MSLP)" begin

nep=nep_gallery("dep0");
nep = DEP(nep.n, nep.A, [0, 0.8])
#nep=shift_and_scale(nep;shift=0,scale=2);

n=size(nep,1);

(λ,v)=newton(nep,v=ones(n))

v=v/norm(v);
S0=reshape([λ],1,1);
V0=reshape(v,n,1);
dnep=effenberger_deflation(nep,S0,V0)

for k=1:6
    λ0=-0.1+0.1im;

    (λ2,v2)=mslp(dnep,maxit=1000,
                        λ=λ0,
                        displaylevel=0
                        );
    println("Computing eigenvalue k=",k);
    v2=v2/norm(v2);
    # Expand the partial Schur factorization with the computed solution
    S0=hcat(S0,v2[(n+1):end]);
    S0=vcat(S0,hcat(zeros(size(S0,1))',λ2))
    V0=hcat(V0,v2[1:n]);
    dnep=effenberger_deflation(nep,S0,V0)
end

# Transform partial Schur factorization to eigenpairs
λv,VS=eigen(S0)
V=V0*VS;

# Test that λv and V are now eigpairs
for i=1:size(λv,1)
    v=V[:,i]/norm(V[:,i]);
    λ=λv[i];
    @test norm(compute_Mlincomb(nep,λ,v))<sqrt(eps())
end

end
