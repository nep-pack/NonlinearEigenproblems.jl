using NonlinearEigenproblems.NEPCore
using Printf

function nleigs_verify_lambdas(nrlambda, nep::NEP, X, lambda, tol = 1e-5)
    @test length(lambda) == nrlambda
    @info "Found $(length(lambda)) lambdas:"
    for i in eachindex(lambda)
        λ = lambda[i]
        nrm = default_errmeasure(nep)(λ, X[:, i])
        @test nrm < tol
        @info "λ[$i] = $λ (norm = $(@sprintf("%.3g", nrm)))"
    end
end
