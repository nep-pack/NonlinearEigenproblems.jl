export @default_compute_rf_inner_solver
# Try to automatically determine an appropriate scalar compute_rf solver.
macro default_compute_rf_inner_solver(nep)
    return (nep isa PEP) ? PolyeigInnerSolver() : ScalarNewtonInnerSolver();
end


"""
    struct ScalarNewtonInnerSolver <: InnerSolver
    function ScalarNewtonInnerSolver(;tol=eps()*100,maxit=80,bad_solution_allowed=true)

This type can be used with  `compute_rf`. It solves the scalar nonlinear problem
using a Rayleig functional. In contrast to the `NewtonInnerSolver`, this
function does not require an SPMF (or precomputations associated with SPMFs).
"""
struct ScalarNewtonInnerSolver <: InnerSolver
    tol::Float64
    maxit::Int
    bad_solution_allowed::Bool
    function ScalarNewtonInnerSolver(;tol=eps()*100,maxit=80,bad_solution_allowed=true)
        return new(tol,maxit,bad_solution_allowed);
    end
end

function compute_rf(T0::Type{T}, nep::NEP, x, inner_solver::ScalarNewtonInnerSolver;
                        y=x, target=zero(T), λ=target,
                        kwargs...) where T
    λ_iter = T(λ);
    Δλ = T(Inf)
    count = 0
    while (abs(Δλ)>inner_solver.tol) & (count<inner_solver.maxit)
        count = count+1
        # compute function value and derivative
        z1 = compute_Mlincomb(nep, λ_iter, reshape(x,size(nep,1),1))
        z2 = compute_Mlincomb(nep, λ_iter, reshape(x,size(nep,1),1),[T(1)],1)

        Δλ = -dot(y,z1)/dot(y,z2);
        λ_iter += Δλ
    end

    if ((count==inner_solver.maxit) && (! (inner_solver.bad_solution_allowed )))
        throw(NoConvergenceException());
    end

    # Return type is a vector of correct type
    λ_star::Array{T,1} = Array{T,1}(undef, 1)
    if (T <: Real) && (typeof(λ_iter) != T) && (imag(λ_iter)/real(λ_iter) < inner_solver.tol)
        # Looking for a real quantity (AND) iterate is not real (AND) complex part is negligible
        λ_star[1] = real(λ_iter) # Truncate to real
    else
        λ_star[1] = λ_iter
    end
    return λ_star
end


"""
    compute_rf(eltype::Type,nep::NEP,x,inner_solver::InnerSolver;
               y=x, target=zero(T), λ=target,TOL=eps(real(T))*1e3,max_iter=10)

Computes the Rayleigh functional of the `nep`, i.e., computes a vector
``Λ`` of values ``λ``
such that ``y^TM(λ)x=0``, using the procedure specified in `inner_solver`.
The default behaviour consists of a scalar valued
Newton-iteration, and the returned vector has only one element.

The given eltype<:Number is the type of the returned vector.

# Example

This uses `iar` to solve the (scalar) nonlinear problem.

```julia-repl
julia> nep=nep_gallery("dep0");
julia> x=ones(size(nep,1));
julia> s=compute_rf(ComplexF64,nep,x,IARInnerSolver())[1] # Take just first element
-1.186623627962043 - 1.5085094961223182im
julia> x'*compute_Mlincomb(nep,s,x)
-8.881784197001252e-16 + 1.0880185641326534e-14im
```
"""
    function compute_rf(T0::Type{T}, nep::NEP, x, inner_solver::InnerSolver;
                        y=x, target=zero(T), λ=target,
                        TOL=eps(real(T))*1e3, max_iter=30,kwargs...) where T


        pnep=create_proj_NEP(nep);
        n=size(nep,1);
        set_projectmatrices!(pnep,reshape(y,n,1),reshape(x,n,1));

        local λv
        try
            (λv,xv)=inner_solve(inner_solver,T0,pnep,σ=target,λv=[λ]);
        catch e
            # Even return the approximation upon failure
            if (e isa NoConvergenceException)
                λv=e.λ
            else
                rethrow(e)
            end
        end
        if (λv isa Vector)
            return λv
        else
            return [λv];
        end
    end




function compute_rf(::Type{T},nep::NEP,x, inner_solver::PolyeigInnerSolver; y=x, target=zero(T), λ=target,
                    TOL=eps(real(T))*1e3,max_iter=10) where T

    a=zeros(T,size(nep.A,1))
    for i=1:size(nep.A,1)
        a[i]=dot(y,nep.A[i]*x);
    end
    m=size(nep.A,1);
    # Build a companion matrix (maybe bigfloat)
    A=zeros(T,m-1,m-1);
    for k=1:m-1
        if (k<m-1)
            A[k,k+1]=1;
        end
        A[end,k]=-a[k]/a[end]
    end
    pp=eigvals(A);
    if (T <: Real) # If we specify real arithmetic, return real eigvals
        pp=real(pp);
    end

    II=sortperm(abs.(pp .- target))
    return pp[II];
end

# Overload a version of eigvals(::Matrix{BigFloat}) for compute_rf
import LinearAlgebra.eigvals
function eigvals(A::Union{Matrix{BigFloat},Matrix{Complex{BigFloat}}})
    T=eltype(A);
    Tfloat= (T <: Complex) ? (ComplexF64) : (Float64)
    Afloat=Matrix{Tfloat}(A);
    λv,X=eigen(Afloat);
    local λv_bigfloat=Vector{Complex{BigFloat}}(λv);
    # Some Rayleigh quotient iteration steps in bigfloat to improve accuracy:
    for j=1:size(λv,1)
        local s=λv_bigfloat[j];
        z=Vector{Complex{BigFloat}}(X[:,j])
        for f=1:3
            z=(A-s*I)\z; normalize!(z);
            s=z'*A*z;
        end
        λv_bigfloat[j]=s
    end
    if (maximum(abs.((λv_bigfloat-λv)./λv_bigfloat)) > eps()*100)
        @warn("Correction iteration of eigvals for bigfloat, made the eigenvalues move considerably")
    end

    return λv_bigfloat
end
