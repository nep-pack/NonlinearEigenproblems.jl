function nleigs_verify_lambdas(nrlambda, nep::NEP, X, lambda, tol = 1e-5)
    @test length(lambda) == nrlambda

    @printf("Found %d lambdas:\n", length(lambda))
    for i in eachindex(lambda)
        λ = lambda[i]
        nrm = default_errmeasure(nep)(λ, X[:, i])
        @test nrm < tol
        @printf("λ[%d] = %s (norm = %.3g)\n", i, λ, nrm)
    end
end
