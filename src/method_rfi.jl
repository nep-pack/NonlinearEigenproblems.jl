using LinearAlgebra
using Random

export rfi
export rfi_b

"""
    rfi(nep,nept,[λ=0,][errmeasure=default_errmeasure,][tol=eps()*100,][maxit=100,][v=randn,][u=randn,][displaylevel=0,][linsolvecreator=default_linsolvecreator,])

This is an implementation of the two-sided Rayleigh functional Iteration. This method requires the transpose of the NEP, which needs to be provided in nept.

# Example:
julia> nep=nep_gallery("dep0");
julia> nept=DEP([nep.A[1]',nep.A[2]'])
julia> λ,v,u=rfi(nep,nept,v=ones(size(nep,1)))
julia> opnorm(compute_Mder(nep,λ)*v) % v is a right eigenvector
5.4672143489065705e-16
julia> opnorm(u'*compute_Mder(nep,λ)) % u is a left eigenvector
4.1447913221215544e-16

# Reference
*  Algorithm 4 in  Schreiber, Nonlinear Eigenvalue Problems: Newton-type Methods and Nonlinear Rayleigh Functionals, PhD thesis, TU Berlin, 2008

"""
rfi(nep::NEP, nept::NEP; kwargs...) = rfi(ComplexF64,nep, nept,;kwargs...)
function rfi(::Type{T},
            nep::NEP,
            nept::NEP;
            errmeasure::Function=default_errmeasure(nep::NEP),
            tol = eps(real(T))*1000,
            maxit=100,
            λ = zero(T),
            v = randn(size(nep,1)),
            u = randn(size(nep,1)),
            linsolvercreator::Function=default_linsolvercreator,
            displaylevel=0) where {T <: Number}

        err = Inf

        #Ensure type coherence
        λ = T(λ)
        v = Vector{T}(v)
        u = Vector{T}(u)
        #Normalize v and u
        normalize!(v)
        normalize!(u)


        try
            for k=1:maxit
                err = errmeasure(λ,u)

                if(err < tol)
                    return λ,u,v
                end

                @ifd(@printf("Iteration: %2d errmeasure:%.18e \n",k, err))

                local linsolver::LinSolver = linsolvercreator(nep,λ)
                local linsolver_t::LinSolver = linsolvercreator(nept,λ)

                #S1: T(λ_k)x_(k+1) = T'(λ_k)u_(k)
                x = lin_solve(linsolver,compute_Mlincomb(nep,λ,u,[T(1)],1),tol = tol)
                u[:] = normalize(x)

                #S2: (T(λ_k)^H)y_(k+1) = (T'(λ_k)^H)v_(k)
                y = lin_solve(linsolver_t,compute_Mlincomb(nept,λ,v,[T(1)],1),tol = tol)
                v[:] = normalize(y)

                λ_vec = compute_rf(nep, u, y=v)
                λ = closest_to(λ_vec,  λ)
            end
        catch e
            isa(e, SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))
            if (errmeasure(λ,u)>tol)
                u[:] = compute_eigvec_from_eigval_lu(nep, λ, default_linsolvercreator)
                normalize!(u)
                v[:] = compute_eigvec_from_eigval_lu(nept, λ, default_linsolvercreator)
                normalize!(v)
            end
            return (λ,u,v)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,u,err,msg))
end

"""
Two-sided Rayleigh functional Iteration(Bordered version), as given as Algorithm 5 in  "Nonlinear Eigenvalue Problems: Newton-type Methods and
Nonlinear Rayleigh Functionals", by Kathrin Schreiber.
"""
rfi_b(nep::NEP, nept::NEP; kwargs...) = rfi_b(ComplexF64,nep, nept,;kwargs...)
function rfi_b(::Type{T},
            nep::NEP,
            nept::NEP;
            errmeasure::Function=default_errmeasure(nep::NEP),
            tol = eps(real(T))*1000,
            maxit=100,
            λ = zero(T),
            v = randn(size(nep,1)),
            u = randn(size(nep,1)),
            linsolvercreator::Function=default_linsolvercreator,
            displaylevel=1) where {T <: Number}

        err = Inf

        #Ensure type coherence
        λ = T(λ)
        v = Vector{T}(v)
        u = Vector{T}(u)
        #Normalize v and u
        normalize!(v)
        normalize!(u)

        try
            for k=1:maxit
                err = errmeasure(λ,u)

                if(err < tol)
                    return λ,u,v
                end

                @ifd(@printf("Iteration: %2d errmeasure:%.18e \n",k, err))

                #Construct C_k
                C = [compute_Mder(nep,λ,0) compute_Mlincomb(nep,λ,u,[T(1)],1);v'*compute_Mder(nep,λ,1) T(0.0)]
                #local linsolver::LinSolver = BackslashLinSolver(C);

                #C[s;μ] = -[T(λ)u;0]
                #l1 = lin_solve(linsolver,-[compute_Mlincomb(nep,λ,u,[1],0);0],tol = tol);
                l1 = C\-[compute_Mlincomb(nep,λ,u,[T(1)],0);0]
                s = l1[1:end-1]
                u[:] = normalize(u+s)

                #C[t;ν] = -[T(λ)'v;0]
                #l2 = lin_solve(linsolver,-[compute_Mlincomb(nept,λ,v,[1],0);0],tol = tol);
                l2 = C\-[compute_Mlincomb(nept,λ,v,[T(1)],0);0]
                t = l2[1:end-1]
                v[:] = normalize(v+t)

                λ_vec = compute_rf(nep, u, y=v)
                λ = closest_to(λ_vec,  λ)
            end
        catch e
            isa(e, SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))
            if (errmeasure(λ,u)>tol)
                u[:] = compute_eigvec_from_eigval_lu(nep, λ, default_linsolvercreator)
                normalize!(u)
                v[:] = compute_eigvec_from_eigval_lu(nept, λ, default_linsolvercreator)
                normalize!(v)
            end
            return (λ,u,v)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,u,err,msg))
end
