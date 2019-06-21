using ..NEPCore
using ..NEPTypes
using ..LinSolvers
export mslp


"""
     mslp([eltype],nep::NEP;[errmeasure,][tol,][maxit,][λ,][v,][logger,][eigsolvertype::Type][armijo_factor=1,][armijo_max])

Runs the method of successive linear problems. The  method requires the solution of a
generalized eigenvalue problem in every iteration. The method used for the eigenvalue
computation is specified in eigsolvertype.
See [`newton`](@ref) for other parameters.


# Example
Create a rational NEP with SPMFs.
```julia-repl
julia> Av=[ones(3,3),eye(3,3),triu(ones(3,3))];
julia> fv=[S-> S, S -> S^2, S::AbstractArray -> inv(Matrix(S)-eye(S)*10)]
julia> nep=SPMF_NEP(Av,fv)
julia> (λ,v)=mslp(nep)
julia> compute_Mlincomb(nep,λ,v)
3-element Array{Complex{Float64},1}:
 -1.38778e-17+1.65715e-18im
 -5.55112e-17+1.30633e-17im
 -4.16334e-17-1.54436e-17im
```

# References
* A. Ruhe, Algorithms for the nonlinear eigenvalue problem, SIAM J. Numer. Anal. 10 (1973) 674-689



"""
mslp(nep::NEP;params...)=mslp(ComplexF64,nep;params...)
function mslp(::Type{T},
                 nep::NEP;
                 errmeasure::ErrmeasureType = DefaultErrmeasure,
                 tol::Real=eps(real(T))*100,
                 maxit::Integer=100,
                 λ::Number=zero(T),
                 logger=0,
                 eigsolvertype::Type=DefaultEigSolver) where {T<:Number}

    @parse_logger_param!(logger)

    # Ensure types λ is of type T
    λ::T = T(λ)

    # Allocate memory for the eigenvector approximation
    v = zeros(T,size(nep,1));
    d = zeros(T,1);


    err=Inf;

    # Init errmeasure
    ermdata=init_errmeasure(errmeasure,nep);

    # Main loop
    for k=1:maxit

        # solve generalized eigenvalue problem
        solver=eigsolvertype(compute_Mder(nep,λ,0),compute_Mder(nep,λ,1));

        # This will throw an error if the eigenvector is not of correct type
        d[:],v[:] = eig_solve(solver,target=0,nev=1);

        # update eigenvalue
        λ += -d[1] # This will throw an error if the eigval update d is not of correct type

        # Normalize
        normalize!(v)

        # Checck for convergence
        err=estimate_error(ermdata,λ,v)
        push_iteration_info!(logger,k,err=err,λ=λ);

        if (err< tol)
            return (λ,v)
        end

    end

    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(λ,v,err,msg))
end
