export jd

using IterativeSolvers

   """
    jd(nep::ProjectableNEP;orthmethod::Type{T_orth} = DGKS,errmeasure::Function = default_errmeasure(nep::NEP),linsolvercreator::Function=default_linsolvercreator,proj_eig_solver::Function = default_proj_eig_solver(nep::ProjectableNEP),Neig = 1,tol = eps(real(T))*100, maxit = 100, λ = zero(T), v0 = randn(size(nep,1)), displaylevel = 0)

Runs the Jacobi-Davidson method with a projected solver in proj_eig_solver.


# Example
```julia-repl
julia> using NonlinearEigenproblems: NEPSolver, NEPCore, Gallery
julia> nep=nep_gallery("real_quadratic");
julia> λ,v=jd(nep,tol=1e-5,maxit=100);
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
8.85629539860136e-12
```

# References
* The Jacobi-Davidson method for nonlinear eigenvalue problems, Betcke.
"""

jd(nep::NEP;params...) = jd(Complex128,nep;params...)
function jd{T, T_orth<:IterativeSolvers.OrthogonalizationMethod}(::Type{T},
                    nep::ProjectableNEP;
                    orthmethod::Type{T_orth} = DGKS,
                    errmeasure::Function = default_errmeasure(nep::NEP),
                    linsolvercreator::Function=default_linsolvercreator,
                    Neig = 1,
                    tol = eps(real(T))*100,
                    maxit = 100,
                    λ = zero(T),
                    v0 = randn(size(nep,1)),
                    displaylevel = 0,
                    inner_solver_method = NEPSolver.DefaultInnerSolver,
                    isHerm = false)

    n = size(nep,1)
    if (maxit > n)
        warn("maxit = ", maxit, " is larger than size of NEP = ", n,". Setting maxit = size(nep,1)")
        maxit = n
    end

    λ_vec::Array{T,1} = Array{T,1}(Neig)
    u_vec::Array{T,2} = zeros(T,n,Neig)
    u::Array{T,1} = Array{T,1}(v0)/norm(v0)
    pk::Array{T,1} = zeros(T,n)

    proj_nep = create_proj_NEP(nep)

    V_memory::Array{T,2} = zeros(T, size(nep,1), maxit+1)
    V_memory[:,1] = u

    if( !isHerm )
        W_memory = zeros(T, size(nep,1), maxit+1)
        W_memory[:,1] = compute_Mlincomb(nep,λ,u);
        if inner_solver_method == NEPSolver.SGIterInnerSolver
            error("Cannot use SGITER if problem is not Hermitian")
        end
    else
        W_memory = view(V_memory, :, :);
    end

    dummy_vector = zeros(T,maxit+1)
    conveig = 0

    for k=1:maxit
        err=errmeasure(λ,u)
        @ifd(print("Iteration: ", k, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
        if convergence_criterion(err, tol, λ, λ_vec, conveig, T)
            conveig += 1
            λ_vec[conveig] = λ
            u_vec[:,conveig] = u
            if (conveig == Neig)
                return (λ_vec,u_vec)
            end
        end

        # Projected matrices
        V = view(V_memory, :, 1:k); W = view(W_memory, :, 1:k); # extact subarrays, memory-CPU efficient
        v = view(V_memory, :, k+1); w = view(W_memory, :, k+1)  # next vector position
        set_projectmatrices!(proj_nep, W, V)


        # find the eigenvalue with smallest absolute value of projected NEP
        λv,sv = inner_solve(inner_solver_method, proj_nep,
                            j = conveig+1, # For SG-iter
                            λv = zeros(T,conveig+1),
                            σ=0,
                            Neig=conveig+1)
        λ,s = jd_eig_sorter(λv, sv, conveig+1)
        s = s/norm(s)

        # the approximate eigenvector
        u[:] = V*s


        # solve for basis extension using comment on top of page 367 to avoid
        # matrix access. The orthogonalization to u comes anyway since u in V
        pk[:] = compute_Mlincomb(nep,λ,u,[1],1)
        linsolver = linsolvercreator(nep,λ)
        v[:] = lin_solve(linsolver, pk, tol=tol) # M(λ)\pk
        orthogonalize_and_normalize!(V, v, view(dummy_vector, 1:k), orthmethod)

        if( !isHerm )
            w[:] = compute_Mlincomb(nep,λ,u);
            orthogonalize_and_normalize!(W, w, view(dummy_vector, 1:k), orthmethod)
        end
    end

    err=errmeasure(λ,u)
    # Check one last time since space was expanded since prevous
    if convergence_criterion(err, tol, λ, λ_vec, conveig, T)
        conveig += 1
        λ_vec[conveig] = λ
        u_vec[:,conveig] = u
        if (conveig == Neig)
            return (λ_vec,u_vec)
        end
    end

    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(Neig)."
    throw(NoConvergenceException(cat(1,λ_vec[1:conveig],λ),cat(2,u_vec[:,1:conveig],u),err,msg))
end


function convergence_criterion(err, tol, λ, λ_vec, conveig, T)
    # Small error and not already found (expception if it is the first)
    # Exclude eigenvalues in a disc of radius of ϵ^(1/4)
    return (err < tol) && (conveig == 0 ||
           all( abs.(λ .- λ_vec[1:conveig])./abs.(λ_vec[1:conveig]) .> sqrt(sqrt(eps(real(T)))) ) )
end

function jd_eig_sorter(λv::Array{T,1}, V, N) where T <: Number
    NN = min(N, length(λv))
    c = sortperm(abs.(λv))
    λ = λv[c[NN]]
    s = V[:,c[NN]]
    return λ, s
end
