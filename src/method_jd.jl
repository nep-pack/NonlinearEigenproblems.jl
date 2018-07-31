export jd

using IterativeSolvers

   """
    function jd([eltype]], nep::ProjectableNEP; [Neig=1], [tol=eps(real(T))*100], [maxit=100], [λ=zero(T)], [orthmethod=DGKS],  [errmeasure=default_errmeasure], [linsolvercreator=default_linsolvercreator], [v0 = randn(size(nep,1))], [displaylevel=0], [inner_solver_method=NEPSolver.DefaultInnerSolver], [projtype=:PetrovGalerkin])
The function computes eigenvalues using Jacobi-Davidson method, which is a projection method.
The projected problems are solved using a solver spcified through the type `inner_solver_method`.
For numerical stability the basis is kept orthogonal, and the method for orthogonalization is specified by `orthmethod`, see the package `IterativeSolvers.jl`.
The function tries to compute `Neig` number of eigenvalues, and throws a `NoConvergenceException` if it cannot.
The value `λ` and the vector `v0` are initial guesses for an eigenpair. `linsolvercreator` is a function which specifies how the linear system is created and solved.
By default the method uses a Petrov-Galerkin framework, with a trial (left) and test (right) space, hence W^H T(λ) V is the projection considered. By specifying  `projtype` to be `:Galerkin` then W=V.


# Example
```julia-repl
julia> using NonlinearEigenproblems: NEPSolver, NEPCore, Gallery
julia> nep=nep_gallery("dep0",50);
julia> λ,v=jd(nep,tol=1e-5,maxit=20);
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
3.4016933647983415e-8
```

# References
* T. Betcke and H. Voss, A Jacobi-Davidson-type projection method for nonlinear eigenvalue problems. Future Gener. Comput. Syst. 20, 3 (2004), pp. 363-372.
* H. Voss, A Jacobi–Davidson method for nonlinear eigenproblems. In: International Conference on Computational Science. Springer, Berlin, Heidelberg, 2004. pp. 34-41.
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
"""

jd(nep::NEP;params...) = jd(Complex128,nep;params...)
function jd(::Type{T},
            nep::ProjectableNEP;
            orthmethod::Type{T_orth} = IterativeSolvers.DGKS,
            errmeasure::Function = default_errmeasure(nep::NEP),
            linsolvercreator::Function=default_linsolvercreator,
            Neig::Int = 1,
            tol = eps(real(T))*100,
            maxit::Int = 100,
            λ = zero(T),
            v0 = randn(size(nep,1)),
            displaylevel = 0,
            inner_solver_method = NEPSolver.DefaultInnerSolver,
            projtype = :PetrovGalerkin)  where {T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod}

    n = size(nep,1)
    if (maxit > n)
        warn("maxit = ", maxit, " is larger than size of NEP = ", n,". Setting maxit = size(nep,1)")
        maxit = n
    end

    λ::T = T(λ)
    λ_vec::Array{T,1} = Array{T,1}(Neig)
    u_vec::Array{T,2} = zeros(T,n,Neig)
    u::Array{T,1} = Array{T,1}(v0); u[:] = u/norm(u);
    pk::Array{T,1} = zeros(T,n)
    proj_nep = create_proj_NEP(nep)
    dummy_vector::Array{T,1} = zeros(T,maxit+1)
    conveig::Int64 = 0

    V_memory::Array{T,2} = zeros(T, size(nep,1), maxit+1)
    V_memory[:,1] = u
    if( projtype == :PetrovGalerkin ) # Petrov-Galerkin uses a left (test) and a right (trial) space
        W_memory = zeros(T, size(nep,1), maxit+1)
        W_memory[:,1] = compute_Mlincomb(nep,λ,u);
        if inner_solver_method == NEPSolver.SGIterInnerSolver
            error("Need to use 'projtype' :Galerkin in order to use SGITER as inner solver.")
        end
    elseif( projtype == :Galerkin ) # Galerkin uses the same trial and test space
        W_memory = view(V_memory, :, :);
    else
        error("Unsupported 'projtype'. The type '", projtype, "' is not supported.")
    end


    # Initial check for convergence
    err = errmeasure(λ,u)
    @ifd(print("Iteration: ", 0, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
    if convergence_criterion(err, tol, λ, λ_vec, conveig, T)
        conveig += 1
        λ_vec[conveig] = λ
        u_vec[:,conveig] = u
        if (conveig == Neig)
            return (λ_vec,u_vec)
        end
    end

    # Main loop
    for k=1:maxit
        # Projected matrices
        V = view(V_memory, :, 1:k); W = view(W_memory, :, 1:k); # extact subarrays, memory-CPU efficient
        v = view(V_memory, :, k+1); w = view(W_memory, :, k+1)  # next vector position
        set_projectmatrices!(proj_nep, W, V)


        # find the eigenvalue with smallest absolute value of projected NEP
        λv,sv = inner_solve(inner_solver_method, T, proj_nep,
                            j = conveig+1, # For SG-iter
                            λv = zeros(T,conveig+1),
                            σ=zero(T),
                            Neig=conveig+1)
        λ,s = jd_eig_sorter(λv, sv, conveig+1)
        s = s/norm(s)

        # the approximate eigenvector
        u[:] = V*s


        # Check for convergence
        err = errmeasure(λ,u)
        @ifd(print("Iteration: ", k, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
        if convergence_criterion(err, tol, λ, λ_vec, conveig, T)
            conveig += 1
            λ_vec[conveig] = λ
            u_vec[:,conveig] = u
            if (conveig == Neig)
                return (λ_vec,u_vec)
            end
        end


        # solve for basis extension using comment on top of page 367 to avoid
        # matrix access. The orthogonalization to u comes anyway since u in V
        pk[:] = compute_Mlincomb(nep,λ,u,[T(1)],1)
        linsolver = linsolvercreator(nep,λ)
        v[:] = lin_solve(linsolver, pk, tol=tol) # M(λ)\pk
        orthogonalize_and_normalize!(V, v, view(dummy_vector, 1:k), orthmethod)

        if( projtype == :PetrovGalerkin )
            w[:] = compute_Mlincomb(nep,λ,u);
            orthogonalize_and_normalize!(W, w, view(dummy_vector, 1:k), orthmethod)
        end
    end



    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(Neig)."
    throw(NoConvergenceException(cat(1,λ_vec[1:conveig],λ),cat(2,u_vec[:,1:conveig],u),err,msg))
end


function convergence_criterion(err, tol, λ, λ_vec, conveig, T)::Bool
    # Small error and not already found (expception if it is the first)
    # Exclude eigenvalues in a disc of radius of ϵ^(1/4)
    return (err < tol) && (conveig == 0 ||
           all( abs.(λ .- λ_vec[1:conveig])./abs.(λ_vec[1:conveig]) .> sqrt(sqrt(eps(real(T)))) ) )
end

function jd_eig_sorter(λv::Array{T,1}, V, N) where T <: Number
    NN = min(N, length(λv))
    c = sortperm(abs.(λv))
    λ::T = λv[c[NN]]
    s::Array{T,1} = V[:,c[NN]]
    return λ, s
end
