export jd
export jd_betcke
export jd_effenberger

using IterativeSolvers


   """
 For Jacobi-Davidson type solvers see ```jd_betcke``` or ```jd_effenberger```
"""
function jd()
    error("The function `jd` is not implemented. See `jd_betcke` or jd_effenberger`.")
end

   """
    function jd_betcke([eltype]], nep::ProjectableNEP; [Neig=1], [tol=eps(real(T))*100], [maxit=100], [λ=zero(T)], [orthmethod=DGKS],  [errmeasure=default_errmeasure], [linsolvercreator=default_linsolvercreator], [v0 = randn(size(nep,1))], [displaylevel=0], [inner_solver_method=NEPSolver.DefaultInnerSolver], [projtype=:PetrovGalerkin], [target=zero(T)])
The function computes eigenvalues using Jacobi-Davidson method, which is a projection method.
The projected problems are solved using a solver spcified through the type `inner_solver_method`.
For numerical stability the basis is kept orthogonal, and the method for orthogonalization is specified by `orthmethod`, see the package `IterativeSolvers.jl`.
The function tries to compute `Neig` number of eigenvalues, and throws a `NoConvergenceException` if it cannot.
The value `λ` and the vector `v0` are initial guesses for an eigenpair. `linsolvercreator` is a function which specifies how the linear system is created and solved.
The `target` is the center around which eiganvlues are computed.
By default the method uses a Petrov-Galerkin framework, with a trial (left) and test (right) space, hence W^H T(λ) V is the projection considered. By specifying  `projtype` to be `:Galerkin` then W=V.


# Example
```julia-repl
julia> using NonlinearEigenproblems: NEPSolver, NEPCore, Gallery
julia> nep=nep_gallery("dep0",50);
julia> λ,v=jd_betcke(nep,tol=1e-5,maxit=20);
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
3.4016933647983415e-8
```

# References
* T. Betcke and H. Voss, A Jacobi-Davidson-type projection method for nonlinear eigenvalue problems. Future Gener. Comput. Syst. 20, 3 (2004), pp. 363-372.
* H. Voss, A Jacobi–Davidson method for nonlinear eigenproblems. In: International Conference on Computational Science. Springer, Berlin, Heidelberg, 2004. pp. 34-41.
See also
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
"""
jd_betcke(nep::NEP;kwargs...) = jd_betcke(Complex128,nep;kwargs...)
function jd_betcke(::Type{T},
                   nep::ProjectableNEP;
                   maxit::Int = 100,
                   Neig::Int = 1,
                   projtype::Symbol = :PetrovGalerkin,
                   inner_solver_method::Type = NEPSolver.DefaultInnerSolver,
                   orthmethod::Type{T_orth} = IterativeSolvers.DGKS,
                   errmeasure::Function = default_errmeasure(nep::NEP),
                   linsolvercreator::Function = default_linsolvercreator,
                   tol::Number = eps(real(T))*100,
                   λ::Number = zero(T),
                   v0::Vector = randn(size(nep,1)),
                   target::Number = zero(T),
                   displaylevel::Int = 0) where {T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod}
    # Initial logical checks
    n = size(nep,1)
    if (maxit > n)
        error("maxit = ", maxit, " is larger than size of NEP = ", n,".")
    end
    if (projtype != :Galerkin) && projtype != :PetrovGalerkin
        error("Only accepted values of 'projtype' are :Galerkin and :PetrovGalerkin.")
    end
    if (projtype != :Galerkin) && (inner_solver_method == NEPSolver.SGIterInnerSolver)
        error("Need to use 'projtype' :Galerkin in order to use SGITER as inner solver.")
    end

    # Allocations and preparations
    λ::T = T(λ)
    target::T = T(target)
    tol::real(T) = real(T)(tol)
    λ_vec::Vector{T} = Vector{T}(Neig)
    u_vec::Matrix{T} = zeros(T,n,Neig)
    u::Vector{T} = Vector{T}(v0); u[:] = u/norm(u);
    conveig = 0

    # Initial check for convergence
    err = errmeasure(λ,u)
    @ifd(print("Iteration: ", 0, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
    if (err < tol) #Frist check, no other eiganvalues can be converged
        conveig += 1
        λ_vec[conveig] = λ
        u_vec[:,conveig] = u
    end
    if (conveig == Neig)
        return (λ_vec,u_vec)
    end

    # More allocations and preparations
    pk::Vector{T} = zeros(T,n)
    s::Vector{T} = zeros(T,n)
    proj_nep = create_proj_NEP(nep)
    dummy_vector::Vector{T} = zeros(T,maxit+1)

    V_memory::Matrix{T} = zeros(T, size(nep,1), maxit+1)
    V_memory[:,1] = u
    if( projtype == :PetrovGalerkin ) # Petrov-Galerkin uses a left (test) and a right (trial) space
        W_memory = zeros(T, size(nep,1), maxit+1)
        W_memory[:,1] = compute_Mlincomb(nep,λ,u);
    else # Galerkin uses the same trial and test space
        W_memory = view(V_memory, :, :);
    end

    # Main loop
    for k = 1:maxit
        # Projection matrices
        V = view(V_memory, :, 1:k); W = view(W_memory, :, 1:k); # extact subarrays, memory-CPU efficient
        v = view(V_memory, :, k+1); w = view(W_memory, :, k+1); # next vector position

        # Project and find the eigenvalue projected NEP
        set_projectmatrices!(proj_nep, W, V)
        λv,sv = inner_solve(inner_solver_method, T, proj_nep,
                            j = conveig+1, # For SG-iter
                            λv = zeros(T,conveig+1),
                            σ=zero(T),
                            Neig=conveig+1)
        λ,s[1:k] = jd_eig_sorter(λv, sv, conveig+1, target)
        s[:] = s/norm(s) #OBS: Hack-ish since s is initilaized with zeros - Perserve type and memory

        # the approximate eigenvector
        u[:] = V*s[1:k]

        # Check for convergence
        err = errmeasure(λ,u)
        @ifd(print("Iteration: ", k, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
        if (err < tol) && (conveig == 0 ||
                all( abs.(λ .- λ_vec[1:conveig])./abs.(λ_vec[1:conveig]) .> sqrt(sqrt(eps(real(T)))) ) )
            conveig += 1
            λ_vec[conveig] = λ
            u_vec[:,conveig] = u
        end
        if (conveig == Neig)
            return (λ_vec,u_vec)
        end

        # Solve for basis extension using comment on top of page 367 of Betcke
        # and Voss, to avoid matrix access. Orthogonalization to u comes anyway
        # since u in V. OBS: Non-standard in JD-literature
        pk[:] = compute_Mlincomb(nep,λ,u,[one(T)],1)
        linsolver = linsolvercreator(nep,λ)
        v[:] = lin_solve(linsolver, pk, tol=tol) # M(λ)\pk
        orthogonalize_and_normalize!(V, v, view(dummy_vector, 1:k), orthmethod)

        if( projtype == :PetrovGalerkin )
            w[:] = compute_Mlincomb(nep,λ,u);
            orthogonalize_and_normalize!(W, w, view(dummy_vector, 1:k), orthmethod)
        end
    end #End main loop

    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(Neig)."
    throw(NoConvergenceException(cat(1,λ_vec[1:conveig],λ),cat(2,u_vec[:,1:conveig],u),err,msg))
end



function jd_eig_sorter(λv::Vector{T}, V, N, target::T) where T <: Number
    NN = min(N, length(λv))
    c = sortperm(abs.(λv-target))
    λ::T = λv[c[NN]]
    s::Vector{T} = V[:,c[NN]]
    return λ, s
end


   """
function jd_effenberger()
# References
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
See also
* T. Betcke and H. Voss, A Jacobi-Davidson-type projection method for nonlinear eigenvalue problems. Future Gener. Comput. Syst. 20, 3 (2004), pp. 363-372.
* H. Voss, A Jacobi–Davidson method for nonlinear eigenproblems. In: International Conference on Computational Science. Springer, Berlin, Heidelberg, 2004. pp. 34-41.
"""
jd_effenberger(nep::NEP;kwargs...) = jd_effenberger(Complex128,nep;kwargs...)
function jd_effenberger(::Type{T},
                        nep::ProjectableNEP;
                        maxit::Int = 100,
                        Neig::Int = 1,
                        inner_solver_method::Type = NEPSolver.DefaultInnerSolver,
                        orthmethod::Type{T_orth} = IterativeSolvers.DGKS,
                        errmeasure::Function = default_errmeasure(nep::NEP),
                        linsolvercreator::Function = default_linsolvercreator,
                        tol::Number = eps(real(T))*100,
                        λ::Number = zero(T),
                        v0::Vector = randn(size(nep,1)),
                        target::Number = zero(T),
                        displaylevel::Int = 0) where {T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod}
    # Initial logical checks
    n = size(nep,1)
    if (maxit > n)
        error("maxit = ", maxit, " is larger than size of NEP = ", n,".")
    end
    if (inner_solver_method == NEPSolver.SGIterInnerSolver)
        error("_Method SGITER not accepted as inner solver since deflated problem not min-max.")
    end

    # Allocations and preparations
    λ::T = T(λ)
    target::T = T(target)
    tol::real(T) = real(T)(tol)
    conveig = 0
    Λ::Matrix{T} = zeros(T,Neig,Neig)
    X::Matrix{T} = zeros(T,n,Neig)

    # Initial check for convergence
    err = errmeasure(λ,u)
    @ifd(print("Iteration: ", 0, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
    if (err < tol) #Frist check, no other eiganvalues can be converged
        conveig += 1
        Λ[conveig,conveig] = λ
        X[:,conveig] = u
    end
    if (conveig == Neig)
        return (Λ,X)
    end


    error("Not implemented yet.")
end

# function jd_project_and_inner_eigsolve!(u, s, proj_nep::Proj_NEP, deflated_nep::DeflatedNEP, V, W, conveig, target, inner_solver_method, k, T)
#     X = deflated_nep.V0
#     Λ = deflated_nep.S0
#
#     n = size(deflated_nep.orgnep,1)
#     nn = size(Λ,1)
#
#     W1 = W[1:n,:]
#     W2 = W[(n+1):(n+nn),:]
#     V1 = V[1:n,:]
#     V2 = V[(n+1):(n+nn),:]
#     set_projectmatrices!(proj_nep, W1, V1)
#
# conveig = nn
#
#     # find the eigenvalue with smallest absolute value of projected NEP
#     λv,sv = inner_solve(inner_solver_method, T, proj_nep,
#                         j = conveig+1, # For SG-iter
#                         λv = zeros(T,conveig+1),
#                         σ=zero(T),
#                         Neig=conveig+1)
#     λ,s[1:k] = jd_eig_sorter(λv, sv, conveig+1, target)
#     s[:] = s/norm(s)
#
#     # the approximate eigenvector
#     u[:] = V*s[1:k]
#     return λ
# end
#
#
# function jd_inner_linear_solver(deflated_nep::DeflatedNEP, λ::T, linsolver, pk::Vector{T}, tol)::Vector{T} where {T<:Number}
#     # If it is a deflated NEP we solve with a Schur complement strategy such that
#     # the user specified solve of M can be used.
#     # (M, U; X^T, 0)(v1;v2) = (y1;y2)
#     # OBS: Assume minimality index = 1
#     # OBS: Forms the Schur complement. Assume that only a few eigenvalues are deflated
#
#     X = deflated_nep.V0
#     Λ = deflated_nep.S0
#
#     n = size(deflated_nep.orgnep,1)
#     nn = size(Λ,1)
#     pk1 = pk[1:n]
#     pk2 = pk[n+1:(n+nn)]
#     if nn == 0 # Corner case, actually empty
#         U = X[:,[]]
#     else
#         U = (compute_MM(deflated_nep,Λ,X) - compute_MM(deflated_nep,λ.*one(Λ),X) ) / (Λ - λ.*one(Λ))
#     end
#
#     # Precompute some reused entities
#     pk1tilde = lin_solve(linsolver, pk1, tol=tol) # pk1tilde = M^{-1}pk1
#     Z::Matrix{T} = zeros(T,n,nn)
#     for i = 1:nn
#         Z[:,i] = lin_solve(linsolver, vec(U[:,i]), tol=tol) # Z = M^{-1}U
#     end
#     S = -X'*Z #Schur complement
#     v2 = S\(pk2 - X'*pk1tilde)
#     v1 = pk1tilde - Z*v2
#     return vcat(v1,v2)
# end
