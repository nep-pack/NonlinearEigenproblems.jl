function nleigs_verify_lambdas(nrlambda, NLEP, X, lambda, tol = 1e-5)
    @test length(lambda) == nrlambda

    @printf("Found %d lambdas:\n", length(lambda))
    for i in eachindex(lambda)
        λ = lambda[i]
        #M = NLEP["B"][1] + λ*NLEP["B"][2] + NLEP["f"][1](λ)*NLEP["C"][1] + NLEP["f"][2](λ)*NLEP["C"][2]
        M = funM(NLEP, λ)
        v = X[:, i]
        nrm = norm(M*v)
        @test nrm < tol
        @printf("λ[%d] = %s (norm = %.3g)\n", i, λ, nrm)
    end
end

function funM(NLEP, λ)
    M = copy(NLEP["B"][1])
    for j = 2:length(NLEP["B"])
        # TODO: M .+= ...
        M += λ^(j-1) * NLEP["B"][j]
    end
    for j = 1:length(NLEP["C"])
        M += NLEP["f"][j](λ) * NLEP["C"][j]
    end
    return M
end
