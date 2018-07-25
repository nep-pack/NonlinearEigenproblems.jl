# Run tests on block Newton

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using IterativeSolvers
    using Base.Test
end

@testset "blocknewton" begin
    nep=nep_gallery("dep0",4);
    V=eye(size(nep,1),3);
    S0=zeros(3,3);
    S,V=blocknewton(nep,S=S0,X=V,displaylevel=1,
                    armijo_factor=0.5,
                    maxit=20)

    λv=eigvals(S);
    for λ=λv
        @test minimum(svdvals(compute_Mder(nep,λ)))<sqrt(eps())
    end

end
