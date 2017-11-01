export sgiter

"""
    λ,v = sgiter([eltype],nep::NEP,j::Integer;[λ_min,][λ_max,][λ,][errmeasure,][tol,][maxit,][displaylevel,][eig_solver,])

Finds the ``j``th eigenvalue of the NEP using self-guarded iteration.
The method only works for Hermitian problems, and the eigenvalues are assumed to be real.
If an interval [λ_min,λ_max] is given, then the Rayleigh functional is assumed to be unique on the interval.
If no interval is given, then the minimum solution is always taken.

"""
sgiter(nep::NEP, j::Integer; params...) = sgiter(Complex128, nep, j; params...)
function sgiter{T}(::Type{T},
                   nep::NEP,
                   j::Integer;
                   λ_min = NaN,
                   λ_max = NaN,
                   λ = zero(real(T)),
                   errmeasure::Function = default_errmeasure(nep),
                   tol = eps(real(T)) * 100,
                   maxit = 100,
                   displaylevel = 0,
                   eigsolvertype::DataType = DefaultEigSolver
                   )

    n = size(nep,1)
    if (j > n) || (j <= 0)
       error("j must be a number between 1 and size(nep) = ", n, ". You have selected j = ", j, ".")
    end
    real_T = real(T)
    λ_min::real_T = real_T(λ_min)
    λ_max::real_T = real_T(λ_max)
    if !( (isnan(λ_min) && isnan(λ_max)) || (!isnan(λ_min) && !isnan(λ_max)) )
        error("A proper interval is not chosen.") #Both should be either given or not given
    end
    if (λ_max < λ_min)
       error("The interval cannot be empty, λ_max >= λ_min is required. Supplied interval [λ_min,λ_max] = [", λ_min, ",", λ_max, "].")
    end
    λ::real_T = real_T(λ)
    if ( (λ < λ_min) || (λ > λ_max) )
       error("The starting guess is outside the interval.")
    end

    v::Array{T,1} = zeros(T,n)
    err = 0

    for k = 1:maxit
       eig_solver = eigsolvertype(compute_Mder(nep, λ, 0))
       v[:] = compute_jth_eigenvector(eig_solver, nep, λ, j)
       λ_vec = compute_rf(real_T, nep, v, TOL = tol/10)
       @ifd(println("compute_rf: ", λ_vec))
       λ = choose_correct_eigenvalue_from_rf(λ_vec, λ_min, λ_max)
       @ifd(println(" λ = ", λ))
       err = errmeasure(λ, v)
       @ifd(print("Iteration:", k, " errmeasure:", err, "\n"))
       if (err < tol)
           return (λ, v)
       end
    end

    msg = "Number of iterations exceeded. maxit = $(maxit)."
    throw(NoConvergenceException(λ, v, err ,msg))
end


function compute_jth_eigenvector(eig_solver, nep::NEP, λ, j)
    n = size(nep,1)
    Λ, V = eig_solve(eig_solver, nev = n)
    p = sortperm(Λ);
    return vec(V[:,p[j]])
end


function choose_correct_eigenvalue_from_rf(λ_vec, λ_min, λ_max)
    if ( isnan(λ_min) && isnan(λ_max) )
        return minimum(λ_vec)
    else
        idxes = (λ_vec .<= λ_max) .& (λ_vec .>= λ_min)
        if sum(idxes) > 1
            error("Multiple values of λ found in the interval.")
        end
        if sum(idxes) == 0
            error("No λ found in the prescribed interval.")
        end
        idx = find(idxes)[1]::Integer #This vector is only 1 element
        return λ_vec[idx]
    end
end
