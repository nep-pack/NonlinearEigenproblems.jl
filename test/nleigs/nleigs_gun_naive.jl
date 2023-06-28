#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test

@bench @testset "NLEIGS: Naive Gun" begin
    nep = nep_gallery("nlevp_native_gun")
    square = [-1-1im; -1+1im; 1+1im; 1-1im]
    λ, X, _ = nleigs(nep, 150^2 .+ 200.0*square, logger=displaylevel)
    verify_lambdas(1, nep, λ, X)
end
