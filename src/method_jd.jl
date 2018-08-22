export jd
export jd_betcke
export jd_effenberger

using IterativeSolvers
using NEPTypes.DeflatedNEP

import NEPTypes.create_proj_NEP
import NEPTypes.set_projectmatrices!
import Base.size
import NEPCore.compute_Mder
#import NEPCore.compute_Mlincomb
#import NEPCore.compute_MM



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
    @ifd(@printf("Iteration: %2d  converged eigenvalues: %2d  errmeasure: %.18e\n", 0, 0, err))
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
        W_memory[:,1] = compute_Mlincomb(nep,λ,u); W_memory[:,1] = W_memory[:,1]/norm(W_memory[:,1])
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
        @ifd(@printf("Iteration: %2d  converged eigenvalues: %2d  errmeasure: %.18e\n", k, conveig, err))
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
    u::Vector{T} = Vector{T}(v0); u[:] = u/norm(u)
    target::T = T(target)
    tol::real(T) = real(T)(tol)
    conveig = 0
    Λ::Matrix{T} = zeros(T,Neig,Neig)
    X::Matrix{T} = zeros(T,n,Neig)

    # Initial check for convergence
    err = errmeasure(λ,u)
    @ifd(@printf("Iteration: %3d  converged eigenvalues: %2d  errmeasure: %.18e\n", 0, 0, err))
    if (err < tol) #Frist check, no other eiganvalues can be converged
        conveig += 1
        Λ[1,1] = λ
        X[:,1] = u
    end

    deflated_nep = effenberger_deflation(nep, Λ[1:conveig,1:conveig], X[:,1:conveig])
    tot_nrof_its = 0

    while true
        # Check for fulfillment. If so, compute the eigenpairs
        if (conveig == Neig)
            λ_vec,uv = eig(Λ)
            u_vec = X * uv
            return λ_vec, u_vec
        end

        λ, u, tot_nrof_its = jd_effenberger_inner(T, deflated_nep, maxit, tot_nrof_its, conveig, inner_solver_method, orthmethod, errmeasure, linsolvercreator, tol, target, displaylevel, Neig)
        conveig += 1 #OBS: minimality index = 1, hence only exapnd by one

        # Expand the partial Schur factorization with the computed solution
        Λ[1:(conveig-1),conveig] = u[(n+1):(n+conveig-1)]
        Λ[conveig,conveig] = λ
        X[:,conveig] = u[1:n]

        deflated_nep = effenberger_deflation(nep, Λ[1:conveig,1:conveig], X[:,1:conveig])
    end

    error("This should not be possible.")
end



function jd_effenberger_inner(::Type{T},
                              deflated_nep::DeflatedNEP,
                              maxit::Int,
                              nrof_its::Int,
                              conveig::Int,
                              inner_solver_method::Type,
                              orthmethod::Type,
                              errmeasure::Function,
                              linsolvercreator::Function,
                              tol::Number,
                              target::Number,
                              displaylevel::Int,
                              Neig::Int) where {T<:Number}
    # Allocations and preparations
    X = deflated_nep.V0
    Λ = deflated_nep.S0
    orgnep = deflated_nep.orgnep

    n = size(orgnep,1) # Size of original problem
    m = size(Λ,1) # Size of deflated subspace

    λ::T = T(target)
    u::Vector{T} = rand(T,n+m); u[:] = u/norm(u)

    pk::Vector{T} = zeros(T,n+m)
    s::Vector{T} = zeros(T,maxit+1-nrof_its)
    proj_nep = create_proj_NEP(deflated_nep)
    dummy_vector::Vector{T} = zeros(T,maxit+1-nrof_its)

    # TODO: This memory allocation can potentially be moved out to outer function for opimization
    V_memory::Matrix{T} = zeros(T, n+m, maxit+1-nrof_its)
    V_memory[:,1] = u
    W_memory::Matrix{T} = zeros(T, n+m, maxit+1-nrof_its)
    W_memory[:,1] = compute_Mlincomb(deflated_nep,λ,u); W_memory[:,1] = W_memory[:,1]/norm(W_memory[:,1])

    err = Inf
    for loop_counter = (nrof_its+1):maxit
        k = loop_counter - nrof_its # Which index are we doing on THIS specific level of deflation
        V = view(V_memory, :, 1:k); W = view(W_memory, :, 1:k); # extact subarrays, memory-CPU efficient
        v = view(V_memory, :, k+1); w = view(W_memory, :, k+1); # next vector position

        # Project and solve the projected NEP
        set_projectmatrices!(proj_nep, W, V)
        λv,sv = inner_solve(inner_solver_method, T, proj_nep,
                            λv = zeros(T,1),
                            σ = zero(T),
                            Neig = 1)
        λ,s[1:k] = jd_eig_sorter(λv, sv, 1, target) #Always closest to target, since deflated
        s[:] = s/norm(s) #OBS: Hack-ish since s is initilaized with zeros - Perserve type and memory

        # the approximate eigenvector
        u[:] = V*s[1:k]

        # Compute residual and check for convergence
        rk = compute_Mlincomb(deflated_nep,λ,u)
        err = norm(rk) # TODO: Can we use a custom error measure?  # TODO: #Error measure, see (9) in Effenberger
        @ifd(@printf("Iteration: %3d  converged eigenvalues: %2d  errmeasure: %.18e  space dimension: %3d\n", loop_counter, conveig, err, k))
        if (err < tol) #Frist check, no other eiganvalues can be converged
            @ifd(print("One eigenvalue converged. Deflating.\n"))
            return (λ, u, loop_counter)
        end

        # Extend the basis
        # The following analogy can be made if w is chosen as [y_k; v_k] in (26) in Effenberger
        # Solve for basis extension using comment on top of page 367 of Betcke
        # and Voss, to avoid matrix access. Orthogonalization to u comes anyway
        # since u in V. OBS: Non-standard in JD-literature
        pk = compute_Mlincomb(deflated_nep,λ,u,[one(T)],1)
        linsolver = linsolvercreator(orgnep, λ)
        jd_inner_effenberger_linear_solver!(v, deflated_nep, λ, linsolver, pk, tol)
        orthogonalize_and_normalize!(V, v, view(dummy_vector, 1:k), orthmethod)
        w[:] = rk
        orthogonalize_and_normalize!(W, w, view(dummy_vector, 1:k), orthmethod)



    end

    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(Neig)."
    # Compute the eigenvalues we have and throw these.
    λ_vec,uv = eig(Λ)
    u_vec = X * uv
    throw(NoConvergenceException(λ_vec, u_vec, err, msg))

end



function compute_U(orgnep, μ, X, Λ)
    if size(Λ,1) == 0 # Corner case, actually empty
        U = X[:,[]]
    else
        μI = μ*one(Λ)
        TXΛ = compute_TXΛ(orgnep, Λ, X)
        U = (TXΛ - compute_MM(orgnep, μI, X))/(Λ - μI)
    end
    return U
end
function compute_Uv(orgnep, μ, X, Λ, v)
    μI = μ*one(Λ)
    t = (Λ - μI)\v
    TXΛ = compute_TXΛ(orgnep, Λ, X)
    return TXΛ*t - compute_MM(orgnep, μI, X)*t
end


compute_TXΛ(deflated_nep::DeflatedNEP, Λ, X) = compute_TXΛ(deflated_nep.orgnep, Λ, X)
function compute_TXΛ(orgnep::NEP, Λ, X)
    return zero(X) # If X and Λ is an ivariant pair, then this block is zero. OBS: Assumed also in the derivation of the algorithm, see Effenberger Lemma 3.1
    # if size(Λ,1) == 0 # Corner case, actually empty
    #     TXΛ = X[:,[]]
    # else
    #     TXΛ = compute_MM(orgnep, Λ, X)
    # end
    # return TXΛ
end


function jd_inner_effenberger_linear_solver!(v, deflated_nep::DeflatedNEP, λ::T, linsolver, pk::Vector{T}, tol) where{T}
    # If it is a deflated NEP we solve with a Schur complement strategy such that
    # the user specified solve of M can be used.
    # (M, U; X^T, 0)(v1;v2) = (y1;y2)
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
    U = compute_U(orgnep, λ, X, Λ)

    # Precompute some reused entities
    pk1tilde = lin_solve(linsolver, pk1, tol=tol) # pk1tilde = M^{-1}pk1
    Z::Matrix{T} = zeros(T,n,m)
    for i = 1:m
        Z[:,i] = lin_solve(linsolver, vec(U[:,i]), tol=tol) # Z = M^{-1}U
    end
    S = -X'*Z #Schur complement
    v2[:] = S\(pk2 - X'*pk1tilde)
    v1[:] = pk1tilde - Z*v2

    # M = compute_Mder(deflated_nep,λ,0)
    # vv = M\pk
    # println(norm(v-vv))
    return v
end



mutable struct JD_Inner_Effenberger_Projected_NEP <: Proj_NEP
    orgnep
    org_proj_nep
    X
    Λ
    V1
    V2
    W1
    W2
    function JD_Inner_Effenberger_Projected_NEP(nep::DeflatedNEP)
        org_proj_nep = create_proj_NEP(nep.orgnep)
        this = new(nep, org_proj_nep, nep.V0, nep.S0)
        return this
    end
end


function create_proj_NEP(deflated_nep::DeflatedNEP)
    return JD_Inner_Effenberger_Projected_NEP(deflated_nep::DeflatedNEP)
end


function set_projectmatrices!(nep::JD_Inner_Effenberger_Projected_NEP, W, V)
    n = size(nep.X,1)
    m = size(nep.Λ,1)
    nep.V1 = V[1:n, :]
    nep.V2 = V[(n+1):(n+m), :]
    nep.W1 = W[1:n, :]
    nep.W2 = W[(n+1):(n+m), :]
    set_projectmatrices!(nep.org_proj_nep, nep.W1, nep.V1)
end


function size(nep::JD_Inner_Effenberger_Projected_NEP,dim=-1)
    n = size(nep.W1, 2);
    if (dim == -1)
        return (n,n)
    else
        return n
    end
end


compute_Mder(nep::JD_Inner_Effenberger_Projected_NEP,λ::Number)=compute_Mder(nep,λ,0)
function compute_Mder(nep::JD_Inner_Effenberger_Projected_NEP,λ::Number,i::Integer)
# THIS IS WRONG!!!
# TODO: Need to compute derivative of U and remove the third term for higher derivative computations
    W1T_M_V1 = compute_Mder(nep.org_proj_nep,λ,i)
    W1T_U_V2 = nep.W1' * compute_U(nep.orgnep.orgnep, λ, nep.X, nep.Λ) * nep.V2
    W2T_A_V1 = nep.W2' * nep.X' * nep.V1 #OBS: Here we assume minimality index = 1
    return  W1T_M_V1 + W1T_U_V2 + W2T_A_V1
end
