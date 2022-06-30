export jd
export jd_betcke
export jd_effenberger

using IterativeSolvers
using LinearAlgebra
using Random
using ..NEPTypes:DeflatedNEP
using ..NEPTypes:deflated_nep_compute_Q




"""
 For Jacobi-Davidson type solvers see
[`jd_betcke`](@ref) or
[`jd_effenberger`](@ref).
"""
function jd()
    error("The function `jd` is not implemented. See `jd_betcke` or jd_effenberger`.")
end


"""
    jd_betcke([eltype]], nep::ProjectableNEP; [neigs=1], [tol=eps(real(T))*100], [maxit=100], [λ=zero(T)], [orthmethod=DGKS],  [errmeasure], [linsolvercreator=DefaultLinSolverCreator()], [v = randn(size(nep,1))], [logger=0], [inner_logger=0], [inner_solver_method=DefaultInnerSolver()], [projtype=:PetrovGalerkin], [target=zero(T)])
The function computes eigenvalues using Jacobi-Davidson method, which is a projection method.
The projected problems are solved using a solver spcified through the type `inner_solver_method`.
The logging of the inner solvers are descided by `inner_logger`, which works in the same way as `logger`.
For numerical stability the basis is kept orthogonal, and the method for orthogonalization is specified by `orthmethod`, see the package `IterativeSolvers.jl`.
The function tries to compute `neigs` number of eigenvalues, and throws a `NoConvergenceException` if it cannot.
The value `λ` and the vector `v` are initial guesses for an eigenpair. `linsolvercreator` is a function which specifies how the linear system is created and solved.
The `target` is the center around which eiganvlues are computed.
By default the method uses a Petrov-Galerkin framework, with a trial (left) and test (right) space, hence ``W^H T(λ) V`` is the projection considered. By specifying  `projtype` to be `:Galerkin` then `W=V`.

See [`augnewton`](@ref) for other parameters.


# Example
```julia-repl
julia> nep=nep_gallery("dep0",50);
julia> λ,v=jd_betcke(nep,tol=1e-8,maxit=20);
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
1.9570222439914163e-10
```

# References
* T. Betcke and H. Voss, A Jacobi-Davidson-type projection method for nonlinear eigenvalue problems. Future Gener. Comput. Syst. 20, 3 (2004), pp. 363-372.
* H. Voss, A Jacobi–Davidson method for nonlinear eigenproblems. In: International Conference on Computational Science. Springer, Berlin, Heidelberg, 2004. pp. 34-41.
See also
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
"""
jd_betcke(nep::NEP;kwargs...) = jd_betcke(ComplexF64,nep;kwargs...)
function jd_betcke(::Type{T},
                   nep::ProjectableNEP;
                   maxit::Int = 100,
                   neigs::Int = 1,
                   projtype::Symbol = :PetrovGalerkin,
                   inner_solver_method = DefaultInnerSolver(),
                   orthmethod = IterativeSolvers.DGKS(),
                   errmeasure::ErrmeasureType = DefaultErrmeasure(nep),
                   linsolvercreator=DefaultLinSolverCreator(),
                   tol::Number = eps(real(T))*100,
                   λ::Number = zero(T),
                   v::Vector = randn(size(nep,1)),
                   target::Number = zero(T),
                   logger = 0,
                   inner_logger = 0) where {T<:Number}

    @parse_logger_param!(logger)
    @parse_logger_param!(inner_logger)

    # Initial logical checks
    n = size(nep,1)
    if (maxit > n)
        error("maxit = ", maxit, " is larger than size of NEP = ", n,".")
    end
    if (projtype != :Galerkin) && projtype != :PetrovGalerkin
        error("Only accepted values of 'projtype' are :Galerkin and :PetrovGalerkin.")
    end
    if (projtype != :Galerkin) && isa(inner_solver_method,SGIterInnerSolver)
        error("Need to use 'projtype' :Galerkin in order to use SGITER as inner solver.")
    end

    # Allocations and preparations
    λ::T = T(λ)
    target::T = T(target)
    tol::real(T) = real(T)(tol)
    λ_vec::Vector{T} = Vector{T}(undef,neigs)
    u_vec::Matrix{T} = zeros(T,n,neigs)
    u::Vector{T} = Vector{T}(v)
    normalize!(u)
    conveig = 0


    # Initial check for convergence
    err = estimate_error(errmeasure,λ,u)
    if (err < tol) #Frist check, no other eiganvalues can be converged
        conveig += 1
        λ_vec[conveig] = λ
        u_vec[:,conveig] = u
    end
    if (conveig == neigs)
        return (λ_vec,u_vec)
    end

    # More allocations and preparations
    pk::Vector{T} = zeros(T,n)
    s_memory::Vector{T} = zeros(T,n)
    proj_nep = create_proj_NEP(nep,maxit,T)
    dummy_vector::Vector{T} = zeros(T,maxit+1)

    V_memory::Matrix{T} = zeros(T, size(nep,1), maxit+1)
    V_memory[:,1] = u
    if( projtype == :PetrovGalerkin ) # Petrov-Galerkin uses a left (test) and a right (trial) space
        W_memory = zeros(T, size(nep,1), maxit+1)
        W_memory[:,1] = compute_Mlincomb(nep,λ,u)
        normalize!(view(W_memory,:,1))
    else # Galerkin uses the same trial and test space
        W_memory = view(V_memory, :, :);
    end

    # Main loop
    for k = 1:maxit
        # Projection matrices
        V = view(V_memory, :, 1:k); W = view(W_memory, :, 1:k); # extact subarrays, memory-CPU efficient
        v = view(V_memory, :, k+1); w = view(W_memory, :, k+1); # next vector position
        s = view(s_memory,1:k)

        # Project and find the eigenvalue projected NEP
        expand_projectmatrices!(proj_nep, W, V)
        λv,sv = inner_solve(inner_solver_method, T, proj_nep,
                            j = conveig+1, # For SG-iter
                            λv = λ .* ones(T,conveig+1),
                            σ = target,
                            neigs = conveig+1,
                            inner_logger = inner_logger)
        λ,s = jd_eig_sorter(λv, sv, conveig+1, target)
        normalize!(s)

        # the approximate eigenvector
        u[:] = V*s

        # Check for convergence
        err = estimate_error(errmeasure,λ,u)
        push_iteration_info!(logger,k,λ=λ,err=err,continues=true)
        push_info!(logger," conveig=$conveig");

        if (err < tol) && (conveig == 0 ||
                all( abs.(λ .- λ_vec[1:conveig])./abs.(λ_vec[1:conveig]) .> sqrt(sqrt(eps(real(T)))) ) )
            conveig += 1
            λ_vec[conveig] = λ
            u_vec[:,conveig] = u
        end
        if (conveig == neigs)
            return (λ_vec,u_vec)
        end

        # Solve for basis extension using comment on top of page 367 of Betcke
        # and Voss, to avoid matrix access. Orthogonalization to u comes anyway
        # since u in V. OBS: Non-standard in JD-literature
        pk[:] = compute_Mlincomb(nep,λ,u,[one(T)],1)
        linsolver::LinSolver=create_linsolver(linsolvercreator,nep,λ)
        v[:] = lin_solve(linsolver, pk, tol=tol) # M(λ)\pk
        orthogonalize_and_normalize!(V, v, view(dummy_vector, 1:k), orthmethod)

        if( projtype == :PetrovGalerkin )
            w[:] = compute_Mlincomb(nep,λ,u);
            orthogonalize_and_normalize!(W, w, view(dummy_vector, 1:k), orthmethod)
        end
    end #End main loop

    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(neigs)."
    throw(NoConvergenceException(cat(λ_vec[1:conveig], λ, dims = 1), cat(u_vec[:,1:conveig], u, dims = 2), err, msg))
end



function jd_eig_sorter(λv::Vector{T}, V, N, target::T) where T <: Number
    NN = min(N, length(λv))
    c = sortperm(abs.(λv.-target))
    λ::T = λv[c[NN]]
    s::Vector{T} = V[:,c[NN]]
    return λ, s
end



"""
    jd_effenberger([eltype]], nep::ProjectableNEP; [maxit=100], [neigs=1], [inner_solver_method=DefaultInnerSolver()], [orthmethod=DGKS()], [linsolvercreator=DefaultLinSolverCreator()], [tol=eps(real(T))*100], [λ=zero(T)], [v = rand(T,size(nep,1))], [target=zero(T)],  [logger=0], [inner_logger=0])
The function computes eigenvalues using the Jacobi-Davidson method, which is a projection method.
Repreated eigenvalues are avoided by using deflation, as presented in the reference by Effenberger.
The projected problems are solved using a solver spcified through the type `inner_solver_method`.
The logging of the inner solvers are descided by `inner_logger`, which works in the same way as `logger`.
For numerical stability the basis is kept orthogonal, and the method for orthogonalization is specified by `orthmethod`, see the package `IterativeSolvers.jl`.
The function tries to compute `neigs` number of eigenvalues, and throws a `NoConvergenceException` if it cannot.
The value `λ` and the vector `v` are initial guesses for an eigenpair. `linsolvercreator` is a function which specifies how the linear system is created and solved.
The `target` is the center around which eiganvalues are computed. For further specifications on the `deflation_mode`, see the function `deflate_eigpair`.

See [`augnewton`](@ref) for other parameters.


# Example
```julia-repl
julia> nep=nep_gallery("dep0",100);
julia> λ,v=jd_effenberger(nep,maxit=30,v=ones(size(nep,1)),λ=0);
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
9.577305772525487e-15
```

# References
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
See also
* T. Betcke and H. Voss, A Jacobi-Davidson-type projection method for nonlinear eigenvalue problems. Future Gener. Comput. Syst. 20, 3 (2004), pp. 363-372.
* H. Voss, A Jacobi–Davidson method for nonlinear eigenproblems. In: International Conference on Computational Science. Springer, Berlin, Heidelberg, 2004. pp. 34-41.
"""
jd_effenberger(nep::NEP;kwargs...) = jd_effenberger(ComplexF64,nep;kwargs...)
function jd_effenberger(::Type{T},
                        nep::ProjectableNEP;
                        maxit::Int = 100,
                        neigs::Int = 1,
                        inner_solver_method = DefaultInnerSolver(),
                        orthmethod = IterativeSolvers.DGKS(),
                        linsolvercreator=DefaultLinSolverCreator(),
                        tol::Number = eps(real(T))*100,
                        λ::Number = rand(T),
                        v::Vector = rand(T,size(nep,1)),
                        target::Number = zero(T),
                        deflation_mode = :Auto,
                        logger = 0,
                        inner_logger = 0) where {T<:Number}

    @parse_logger_param!(logger)
    @parse_logger_param!(inner_logger)

    # Initial logical checks
    n = size(nep,1)
    if (maxit > n)
        error("maxit = ", maxit, " is larger than size of NEP = ", n,".")
    end
    if isa(inner_solver_method,SGIterInnerSolver)
        error("Inner solver 'SGIterInnerSolver' not accepted since deflated problem not min-max.")
    end

    # Allocations and preparations
    λ::T = T(λ)
    u::Vector{T} = Vector{T}(v)
    normalize!(u)
    λ_init::T = λ
    u_init::Vector{T} = u

    target::T = T(target)
    tol::real(T) = real(T)(tol)
    conveig = 0
    tot_nrof_its = 0
    V_memory_base::Matrix{T} = zeros(T, n+neigs, maxit+1)
    W_memory_base::Matrix{T} = zeros(T, n+neigs, maxit+1)

    # Initial check for convergence
    err = norm(compute_Mlincomb(nep,λ,u))
    if (err < tol) # Check if initial guess is good enough, otherwise compute initial invariant pair
        λ_init = T(rand())
        u_init = rand(n+1)
    else
        V_memory = view(V_memory_base, 1:n, tot_nrof_its+1:(maxit+1))
        W_memory = view(W_memory_base, 1:n, tot_nrof_its+1:(maxit+1))
        λ, u, tot_nrof_its, u_init, λ_init = dispatch_inner!(T, V_memory, W_memory,
                                                          nep, maxit, tot_nrof_its,
                                                          conveig, inner_solver_method, orthmethod,
                                                          linsolvercreator, tol, target, logger,
                                                          neigs, u_init, λ_init, inner_logger)
    end

    conveig += 1
    deflated_nep = deflate_eigpair(nep, λ, u, mode=deflation_mode)


    while true # Can only escape the loop on convergence (return) or too many iterations (error)
        # Check for fulfillment. If so, compute the eigenpairs and return
        if (conveig == neigs)
            (λ_vec, u_vec,) = get_deflated_eigpairs(deflated_nep)
            return λ_vec, u_vec
        end

        V_memory = view(V_memory_base, 1:(n+conveig), tot_nrof_its+1:(maxit+1))
        W_memory = view(W_memory_base, 1:(n+conveig), tot_nrof_its+1:(maxit+1))
        λ, u, tot_nrof_its, u_init, λ_init = dispatch_inner!(T, V_memory, W_memory,
                                                          deflated_nep, maxit, tot_nrof_its,
                                                          conveig, inner_solver_method, orthmethod,
                                                          linsolvercreator, tol, target, logger,
                                                          neigs, u_init, λ_init, inner_logger)
        conveig += 1 #OBS: minimality index = 1, hence only exapnd by one

        deflated_nep = deflate_eigpair(deflated_nep, λ, u)
    end
end


# Two dispatch functions depending on if it is a DeflatedNEP or a ProjectableNEP
# A bit hacky, but avoids corner cases with "empty deflation"
function dispatch_inner!(::Type{T},
                              V_memory::SubArray,
                              W_memory::SubArray,
                              target_nep::DeflatedNEP,
                              args...) where {T<:Number}
    return  jd_effenberger_inner!(T, V_memory, W_memory,
                                  target_nep, target_nep.V0, target_nep.S0, target_nep.orgnep,
                                  args...)
end
function dispatch_inner!(::Type{T},
                              V_memory::SubArray,
                              W_memory::SubArray,
                              target_nep::ProjectableNEP,
                              args...) where {T<:Number}
    return  jd_effenberger_inner!(T, V_memory, W_memory,
                                  target_nep, zeros(T,size(target_nep,1),0), zeros(T,0,0), target_nep,
                                  args...)
end


function jd_effenberger_inner!(::Type{T},
                              V_memory::SubArray,
                              W_memory::SubArray,
                              target_nep::Union{DeflatedNEP,ProjectableNEP}, # OBS: target_nep can either be a DeflatedNEP, or another NEP if we have no initial invariant pair
                              X::Matrix{T}, # This is explicitly allocated by the dispatch_inner!(...)
                              Λ::Matrix{T}, # This is explicitly allocated by the dispatch_inner!(...)
                              orgnep::NEP, # This is explicitly allocated by the dispatch_inner!(...)
                              maxit::Int,
                              nrof_its::Int,
                              conveig::Int,
                              inner_solver_method,
                              orthmethod::Type,
                              linsolvercreator,
                              tol::Number,
                              target::Number,
                              logger::Logger,
                              neigs::Int,
                              u::Vector{T},
                              λ::T,
                              inner_logger::Logger) where {T<:Number}
    # Allocations and preparations
    n = size(orgnep,1) # Size of original problem
    m = size(Λ,1) # Size of deflated subspace

    normalize!(u)

    λ_temp::T = λ
    newton_step::Vector{T} = Vector{T}(rand(n+m))
    pk::Vector{T} = zeros(T,n+m)
    s_memory::Vector{T} = zeros(T,maxit+1-nrof_its)
    proj_nep = create_proj_NEP(target_nep,maxit+1-nrof_its,T)

    dummy_vector::Vector{T} = zeros(T,maxit+1-nrof_its)

    V_memory[:,1] = u
    W_memory[:,1] = compute_Mlincomb(target_nep, λ, u)
    normalize!(view(W_memory,:,1))

    err = Inf
    for loop_counter = (nrof_its+1):maxit
        k = loop_counter - nrof_its # Which index are we doing on THIS specific level of deflation
        V = view(V_memory, :, 1:k); W = view(W_memory, :, 1:k); # extact subarrays, memory-CPU efficient
        v = view(V_memory, :, k+1); w = view(W_memory, :, k+1); # next vector position
        s = view(s_memory,1:k)

        # Project and solve the projected NEP
        expand_projectmatrices!(proj_nep, W, V)
        λv,sv = inner_solve(inner_solver_method, T, proj_nep,
                            tol = tol/10,
                            λv = λ .* ones(T,2),
                            σ = target,
                            neigs = 2,
                            inner_logger = inner_logger)
        λ_temp,s = jd_eig_sorter(λv, sv, 1, target) #Always closest to target, since deflated
        normalize!(s)

        # If inner solver converged to a solution to the projected problem use that.
        # Otherwise take the "newton step". Not implemented exactly as in Effenbergerm but similar
        if !isnan(λ_temp) && !any(isnan.(s[1:k])) && (norm(compute_Mlincomb(proj_nep, λ_temp, s[1:k])) < tol*50) # Has converged to an eigenvalue of the projected problem
            u[:] = V*s # The approximate eigenvector
            λ = λ_temp
        else
            u[:] = u[:] + newton_step
            normalize!(u)
        end

        # Compute residual and check for convergence
        rk::AbstractVector = compute_Mlincomb(target_nep, λ, u)
        err = norm(rk) #Error measure, see (9)/(3.2) in Effenberger # TODO: Can we use another error measure?
        push_iteration_info!(logger,loop_counter,λ=λ,err=err,continues=true);
        push_info!(logger," conveig=$conveig, subspace dim=$k");

        if (err < tol) #Frist check, no other eiganvalues can be converged
            info_str = "One eigenvalue converged."
            if conveig+1 < neigs
                info_str = string(info_str, " Deflating and restarting.")
            end
            push_info!(logger, info_str)

            # TODO: Here one can implement a continuation with the same basis as in Effenberger section 4.2.5
            # What is here is only a light kind adapted to an "unknown" inner solver
            λ2,s2 = jd_eig_sorter(λv, sv, 2, target)
            if (size(sv,2) > 1) && (abs(λ-λ2)/abs(λ) > sqrt(eps(real(T))))
                normalize!(s2)
                u2 = vcat(V*s2, zero(T))
            else
                λ2 = T(rand())
                u2 = Vector{T}(rand(n+m+1))
            end
            return (λ, u, loop_counter, u2, λ2)

        end

        # Extend the basis
        # The following analogy can be made if w is chosen as [y_k; v_k] in (26)/(4.6) in Effenberger
        # Solve for basis extension using comment on top of page 367 of Betcke
        # and Voss, to avoid matrix access. Orthogonalization to u comes anyway
        # since u in V. OBS: Non-standard in JD-literature
        pk[:] = compute_Mlincomb(target_nep, λ, u, [one(T)], 1)
        linsolver::LinSolver=create_linsolver(linsolvercreator,orgnep,λ)
        jd_inner_effenberger_linear_solver!(v, target_nep, λ, linsolver, pk, tol)
        newton_step = copy(v)
        orthogonalize_and_normalize!(V, v, view(dummy_vector, 1:k), orthmethod)
        w[:] = rk
        orthogonalize_and_normalize!(W, w, view(dummy_vector, 1:k), orthmethod)



    end

    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(neigs)."
    # Compute the eigenvalues we have and throw these.
    F = eigen(Λ)
    λ_vec = F.values
    uv = F.vectors
    u_vec = X * uv
    throw(NoConvergenceException([λ_vec; λ], [u_vec u[1:n]], err, msg))

end


function jd_inner_effenberger_linear_solver!(v, deflated_nep::DeflatedNEP, λ::T, linsolver, pk::Vector{T}, tol) where{T}
    # If it is a deflated NEP we solve with a Schur complement strategy such that
    # the user specified solve of M can be used.
    # (M, U; X^T, 0)(v1;v2) = (pk1;pk2)
    # OBS: Assume minimality index = 1
    # OBS: Forms the Schur complement. Assume that only a few eigenvalues are deflated

    X = deflated_nep.V0
    Λ = deflated_nep.S0
    orgnep = deflated_nep.orgnep
    n = size(orgnep,1)
    m = size(Λ,1)

    v1 = view(v, 1:n)
    v2 = view(v, (n+1):(n+m))
    pk1 = pk[1:n]
    pk2 = pk[(n+1):(n+m)]
    U = deflated_nep_compute_Q(deflated_nep, λ, 0)

    # Precompute some reused entities
    pk1tilde = lin_solve(linsolver, pk1, tol=tol) # pk1tilde = M^{-1}pk1
    Z::Matrix{T} = zeros(T,n,m)
    for i = 1:m
        Z[:,i] = lin_solve(linsolver, vec(U[:,i]), tol=tol) # Z = M^{-1}U
    end
    S = -X'*Z #Schur complement
    v2[:] = S\(pk2 - X'*pk1tilde)
    v1[:] = pk1tilde - Z*v2

    return v
end

function jd_inner_effenberger_linear_solver!(v, nep::NEP, λ::T, linsolver, pk::Vector{T}, tol) where{T}
    v[:] = lin_solve(linsolver, pk, tol=tol)
    return v
end
