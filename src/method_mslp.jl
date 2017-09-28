using NEPCore
using NEPTypes
using LinSolvers
export mslp


"""
     mslp(nep,..)
   Method of successive linear problems
"""
mslp(nep::NEP;params...)=mslp(Complex128,nep;params...)
function mslp{T}(::Type{T},
                 nep::NEP;
                 errmeasure::Function =
                 default_errmeasure(nep::NEP),
                 tol=eps(real(T))*100,
                 maxit=100,
                 λ=zero(T),
                 v=randn(nep.n),
                 displaylevel=0,
                 eigsolvertype::DataType=DefaultEigSolver)

    # Ensure types λ and v are of type T
    λ = T(λ)
    v = Array{T,1}(v)

    σ=λ;     
    err=Inf;


    #levsolver = LinEigSolver();
        
    # Main loop
    for k=1:maxit
        # Normalize
        v=v/norm(v);

        err=errmeasure(λ,v)

        if (displaylevel>0)
            println("Iteration:",k," errmeasure:",err)
        end
        if (err< tol)
            return (λ,v)
        end

        # solve generalized eigenvalue problem
        #d,v = levsolver.solve(compute_Mder(nep,λ,0),B=compute_Mder(nep,λ,1),λ_t=λ,nev=1);
        solver=eigsolvertype(compute_Mder(nep,λ,0),compute_Mder(nep,λ,1));

        d,v = eig_solve(solver,target=λ,nev=1);
        # update eigenvalue
        λ += -d
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
