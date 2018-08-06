function nleigs_verify_lambdas(nrlambda, NLEP, X, lambda, tol = 1e-5)
    @test length(lambda) == nrlambda

    @printf("Found %d lambdas:\n", length(lambda))
    for i in eachindex(lambda)
        λ = lambda[i]
        M = funM(NLEP, λ)
        v = X[:, i]
        nrm = norm(M*v)
        @test nrm < tol
        @printf("λ[%d] = %s (norm = %.3g)\n", i, λ, nrm)
    end
end

function nleigs_verify_lambdas(nrlambda, nep::SPMF_NEP, X, lambda, tol = 1e-5)
    @test length(lambda) == nrlambda

    @printf("Found %d lambdas:\n", length(lambda))
    for i in eachindex(lambda)
        λ = lambda[i]
        nrm = norm(compute_Mlincomb(nep, λ, X[:, i]))
        @test nrm < tol
        @printf("λ[%d] = %s (norm = %.3g)\n", i, λ, nrm)
    end
end

function funM(NLEP, λ)
    M = copy(NLEP["B"][1])
    for j = 2:length(NLEP["B"])
        M += λ^(j-1) * NLEP["B"][j]
    end
    as_matrix(x::Number) = (M = Matrix{eltype(x)}(1,1); M[1] = x; M)
    for j = 1:length(NLEP["C"])
        M += NLEP["f"][j](as_matrix(λ))[1] * NLEP["C"][j]
    end
    return M
end
