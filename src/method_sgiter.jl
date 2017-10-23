export sgiter

"""
    λ,v = sgiter([eltype],nep::NEP,j::Integer;[errmeasure,][tol,][maxit,][displaylevel,][eig_solver,])

Finds the ``j``th eigenvalue of the NEP using self-guarded iteration.

"""
sgiter(nep::NEP, j::Integer; params...) = sgiter(Complex128, nep, j; params...)
function sgiter{T}(::Type{T},
                   nep::NEP,
                   j::Integer;
                   λ = zero(T),
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
    λ::real(T) = real(T(λ))
    local λv::Array{real(T),1};
    v::Array{T,1} = zeros(T,n)
    err = 0

    for k = 1:maxit
       eig_solver = eigsolvertype(compute_Mder(nep, λ, 0))
       println(typeof(eig_solver))
       v[:] = compute_jth_eigenvector(eig_solver, nep, λ, j)
       λv = compute_rf(T, nep, v, TOL = tol/10, max_iter = 10)
       λ=maximum(λv);
       @ifd(println("compute_rf:",λv, " λ=",λ))
       err = errmeasure(λ, v)
       @ifd(print("Iteration:", k, " errmeasure:", err, "\n"))
       if (err < tol)
           return (λ, v)
       end
    end

    msg = "Number of iterations exceeded. maxit = $(maxit)."
    throw(NoConvergenceException(λ, v, err ,msg))
end


function compute_jth_eigenvector(eig_solver::Union{DefaultEigSolver,NativeEigSolver}, nep::NEP, λ, j)
    n = size(nep,1)
    Λ, V = eig_solve(eig_solver, nev = n)
    p = sortperm(Λ);
    return vec(V[:,p[j]])
end
