using NonlinearEigenproblems
using Test
using LinearAlgebra
using IterativeSolvers
@testset "AAAeigs" begin
    nep=nep_gallery("dep0");
    circ=exp.(1im*(0:0.01:2)*pi); circ=circ[1:end-1];
    Z=2*circ;

    v0=normalize(ones(size(nep,1))); # No random

    # Test general nonlinear problem
    (λ,v)=AAAeigs(nep,Z,v0=v0);
    @test length(λ) == 6 # This problem has 6 eigvals smaller than 2.0
    @test count(map(i -> norm(compute_Mlincomb(nep,λ[i],normalize(v[:,i]))),1:6) .< sqrt(eps())) == 6

    # Test SumNEP
    Av=get_Av(nep);
    nep2=SumNEP(PEP([Av[2],-Av[1]]),SPMF_NEP([Av[3]],[s->exp(-s)]));
    (λ,v)=AAAeigs(nep2,Z,weighted=true,v0=v0, return_details=true);
    @test length(λ) == 6
    @test count(map(i -> norm(compute_Mlincomb(nep,λ[i],normalize(v[:,i]))),1:6) .< sqrt(eps())) == 6


end
