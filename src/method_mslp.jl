using ..NEPCore
using ..NEPTypes
using ..LinSolvers
export mslp


"""
     mslp([eltype],nep::NEP;[errmeasure,][tol,][maxit,][λ,][v,][displaylevel,][eigsolvertype::DataType][armijo_factor=1,][armijo_max])

Runs the method of successive linear problems. The  method requires the solution of a
generalized eigenvalue problem in every iteration. The method used for the eigenvalue
computation is specified in eigsolvertype. See `newton` for other parameters.

# Example
Create a rational NEP with SPMFs.
```julia-repl
julia> Av=[ones(3,3),eye(3,3),triu(ones(3,3))];
julia> fv=[S-> S, S -> S^2, S::AbstractArray -> inv(full(S)-eye(S)*10)]
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
mslp(nep::NEP;params...)=mslp(Complex128,nep;params...)
function mslp{T}(::Type{T},
                 nep::NEP;
                 errmeasure::Function =
                 default_errmeasure(nep::NEP),
                 tol=eps(real(T))*100,
                 maxit=100,
                 λ=zero(T),
                 displaylevel=0,
                 eigsolvertype::DataType=DefaultEigSolver)

    # Ensure types λ is of type T
    λ = T(λ)

    # Allocate memory for the eigenvector approximation
    v = zeros(T,size(nep,1));

    err=Inf;

    # Main loop
    for k=1:maxit

        # solve generalized eigenvalue problem
        solver=eigsolvertype(compute_Mder(nep,λ,0),compute_Mder(nep,λ,1));

        d,v = eig_solve(solver,target=0,nev=1);

        # update eigenvalue
        λ += -d

        # Normalize
        v=v/norm(v);

        # Checck for convergence
        err=errmeasure(λ,v)
        @ifd(println("Iteration:",k," errmeasure:",err, " λ=",λ))
        
        if (err< tol)
            return (λ,v)
        end
        
    end

    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(λ,v,err,msg))
end

# A naive implementation of inverse iteration of generalized
# linear eigenvalue problems
function inv_it(nep::NEP,λ=0,v0=randn(nep.n),iters=2)
    Mp=compute_Mder(nep,λ,1);
    M=compute_Mder(nep,λ)
    local w=copy(v0)
    for i=1:iters
        w=M\(Mp*w)
        w=w/norm(w)
    end
    Δ=dot(w,M*w)/dot(w,Mp*w) # Comp delta with Rayleigh Quotient
    return Δ,w
end
