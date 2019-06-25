export sgiter

"""
    λ,v = sgiter([eltype],nep::NEP,j::Integer;[λ_min,][λ_max,][λ,][errmeasure,][tol,][maxit,][logger,][eigsolvertype::Type,])

Finds the `j`-th eigenvalue of the NEP using safeguarded iteration, with eigenvalue numbering
according to min-max theory. The method only works for Hermitian problems, and the eigenvalues
are assumed to be real. If an interval [`λ_min`,`λ_max`] is given, then the Rayleigh functional
is assumed to be unique on the interval. If no interval is given, then the minimum solution is
always taken. The method requires the computation of (all) eigenvalues of a matrix. The `eigsolvertype`
is a `Type` that specifies which eigevalue solver is used inside the algorithm.

See [`newton`](@ref) for other parameters.


# Example
```julia-repl
julia> nep = nep_gallery("real_quadratic");
julia> λ,v = sgiter(nep, 1, λ_min = -10, λ_max = 0,  λ = -10, maxit = 100);
julia> minimum(svdvals(compute_Mder(nep,λ)))
0.0
julia> norm(v)
1.0
```

# References
* V. Mehrmann and H. Voss, Nonlinear eigenvalue problems: a challenge for modern eigenvalue methods, GAMM‐Mitteilungen 27.2 (2004): 121-152.
* H. Voss and B. Werner, Solving sparse nonlinear eigenvalue problems. Technical Report 82/4, Inst. f. Angew. Mathematik, Universität Hamburg, 1982.
* B. Werner. Das Spektrum von Operatorenscharen mit verallgemeinerten Rayleighquotienten. PhD thesis, Fachbereich Mathematik, Universität Hamburg, 1970

"""
sgiter(nep::NEP, j::Integer; params...) = sgiter(ComplexF64, nep, j; params...)
function sgiter(::Type{T},
                   nep::NEP,
                   j::Integer;
                   λ_min::Real = NaN,
                   λ_max::Real = NaN,
                   λ::Number = zero(real(T)),
                   errmeasure::ErrmeasureType = DefaultErrmeasure,
                   tol::Real = eps(real(T)) * 100,
                   maxit::Integer = 100,
                   logger = 0,
                   eigsolvertype::Type = DefaultEigSolver
                   ) where {T<:Number}

    @parse_logger_param!(logger)

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

    # Init errmeasure
    ermdata=init_errmeasure(errmeasure,nep);

    for k = 1:maxit
       eig_solver = eigsolvertype(compute_Mder(nep, λ, 0))
       v[:] = compute_jth_eigenvector(eig_solver, nep, λ, j)
       λ_vec = compute_rf(real_T, nep, v, TOL = tol/10)
       push_info!(logger,2,"compute_rf: $λ_vec")
       λ = choose_correct_eigenvalue_from_rf(λ_vec, λ_min, λ_max)
       err = estimate_error(ermdata,λ, v)
       push_iteration_info!(logger,k,err=err,λ=λ);
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
    if (n>1)
        p = sortperm(Λ)
        return vec(V[:,p[j]])
    else
        return V
    end
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
        idx = findall(idxes)[1]::Int #This vector is only 1 element
        return λ_vec[idx]
    end
end
