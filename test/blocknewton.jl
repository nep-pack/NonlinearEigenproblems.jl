# Run tests on block Newton

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems.NEPCore
using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using Test
using LinearAlgebra

@bench @testset "blocknewton" begin
    nep=nep_gallery("dep0",4);
    V = Matrix(1.0I, size(nep,1), 3)
    S0=zeros(3,3);
    S,V=blocknewton(nep,S=S0,X=V,displaylevel=1,
                    armijo_factor=0.5,
                    maxit=20)

    λv=eigvals(S);
    for λ=λv
        @test minimum(svdvals(compute_Mder(nep,λ)))<sqrt(eps())
    end

end
