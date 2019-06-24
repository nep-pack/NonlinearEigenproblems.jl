using LinearAlgebra
using Random

export rfi
export rfi_b

"""
    rfi(nep,nept,[λ=0,][errmeasure,][tol=eps()*100,][maxit=100,][v=randn,][u=randn,][logger=0,][linsolvecreator=default_linsolvecreator,])


This is an implementation of the two-sided Rayleigh functional Iteration (RFI) to compute an eigentriplet of the problem specified by `nep`.
This method requires the transpose of the NEP, specified in `nept`.
`λ`, `u` and `v` are initial guesses for the eigenvalue, the right eigenvector and the left eigenvector respectively.
See [`newton`](@ref) for other parameters.

# Example
```julia-repl
julia> nep=nep_gallery("dep0");
julia> nept=DEP([nep.A[1]',nep.A[2]'])
julia> λ,v,u=rfi_b(nep,nept)
julia> compute_resnorm(nep,λ,v) % v is a right eigenvector
4.347204570675246e-16
julia> compute_resnorm(nept,λ,u) % u is a left eigenvector
7.173081573164097e-16
```
# Reference
*  Algorithm 4 in  Schreiber, Nonlinear Eigenvalue Problems: Newton-type Methods and Nonlinear Rayleigh Functionals, PhD thesis, TU Berlin, 2008.
"""
rfi(nep::NEP, nept::NEP; kwargs...) = rfi(ComplexF64,nep, nept,;kwargs...)
function rfi(::Type{T},
            nep::NEP,
            nept::NEP;
            errmeasure::ErrmeasureType = DefaultErrmeasure,
            tol = eps(real(T))*1000,
            maxit=100,
            λ::Number = zero(T),
            v::Vector = randn(T,size(nep,1)),
            u::Vector = randn(T,size(nep,1)),
            linsolvercreator::Function=default_linsolvercreator,
            logger=0) where {T <: Number}

        @parse_logger_param!(logger)

        err = Inf

        #Normalize v and u
        normalize!(v)
        normalize!(u)

        # Init errmeasure
        ermdata=init_errmeasure(errmeasure,nep);

        for k=1:maxit
            err = estimate_error(ermdata,λ,u)

            if(err < tol)
                return λ,u,v
            end

            push_iteration_info!(logger,k,err=err,λ=λ,v=v, continues=true);
            push_info!(logger," u=$u");

            local linsolver::LinSolver = linsolvercreator(nep,λ)
            local linsolver_t::LinSolver = linsolvercreator(nept,λ)

            x = lin_solve(linsolver,compute_Mlincomb(nep,λ,u,[T(1)],1),tol = tol)
            u[:] = normalize(x)

            y = lin_solve(linsolver_t,compute_Mlincomb(nept,λ,v,[T(1)],1),tol = tol)
            v[:] = normalize(y)

            λ_vec = compute_rf(nep, u, y=v)
            λ = closest_to(λ_vec,  λ)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,u,err,msg))
end

"""
    rfi_b(nep,nept,[λ=0,][errmeasure,][tol=eps()*100,][maxit=100,][v=randn,][u=randn,][logger=0,][linsolvecreator=default_linsolvecreator,])

This is an implementation of the two-sided Rayleigh functional Iteration(RFI)-Bordered version to compute an eigentriplet of the problem specified by `nep`.
This method requires the transpose of the NEP, specified in `nept`.
`λ`, `u` and `v` are initial guesses for the eigenvalue, the right eigenvector and the left eigenvector respectively.
See [`newton`](@ref) for other parameters.

# Example
```julia-repl
julia> nep=nep_gallery("dep0");
julia> nept=DEP([nep.A[1]',nep.A[2]'])
julia> λ,v,u=rfi_b(nep,nept,v=ones(size(nep,1)))
julia> compute_resnorm(nep,λ,v) % v is a right eigenvector
5.343670589284583e-15
julia> compute_resnorm(nept,λ,u) % u is a left eigenvector
5.271390516634306e-16
```

# Reference
*  Algorithm 5 in  Schreiber, Nonlinear Eigenvalue Problems: Newton-type Methods and Nonlinear Rayleigh Functionals, PhD thesis, TU Berlin, 2008.

"""
rfi_b(nep::NEP, nept::NEP; kwargs...) = rfi_b(ComplexF64,nep, nept,;kwargs...)
function rfi_b(::Type{T},
            nep::NEP,
            nept::NEP;
            errmeasure::ErrmeasureType = DefaultErrmeasure,
            tol = eps(real(T))*1000,
            maxit=100,
            λ::Number = zero(T),
            v::Vector = randn(T,size(nep,1)),
            u::Vector = randn(T,size(nep,1)),
            linsolvercreator::Function=default_linsolvercreator,
            logger=0) where {T <: Number}

        @parse_logger_param!(logger)

        err = Inf
        #Normalize v and u
        normalize!(v)
        normalize!(u)

        # Init errmeasure
        ermdata=init_errmeasure(errmeasure,nep);

        for k=1:maxit
            err = estimate_error(ermdata,λ,u)

            if(err < tol)
                return λ,u,v
            end

            push_iteration_info!(logger,k,err=err,λ=λ,v=v, continues=true);
            push_info!(logger," u=$u");

            #Construct C_k
            C = [compute_Mder(nep,λ,0) compute_Mlincomb(nep,λ,u,[T(1)],1);v'*compute_Mder(nep,λ,1) T(0.0)]

            l1 = C\-[compute_Mlincomb(nep,λ,u,[T(1)],0);0]
            s = l1[1:end-1]
            u[:] = normalize(u+s)

            l2 = C\-[compute_Mlincomb(nept,λ,v,[T(1)],0);0]
            t = l2[1:end-1]
            v[:] = normalize(v+t)

            λ_vec = compute_rf(nep, u, y=v)
            λ = closest_to(λ_vec,  λ)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,u,err,msg))
end
