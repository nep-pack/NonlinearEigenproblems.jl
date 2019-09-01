  """
    compute_rf([eltype],nep::NEP,x; y=x, target=zero(T), λ0=target,TOL=eps(real(T))*1e3,max_iter=10)

Computes the Rayleigh functional of nep, i.e., computes a vector ``Λ`` of values ``λ``
such that ``y^TM(λ)x=0``. The default behaviour consists of a scalar valued
Newton-iteration, and the returned vector has only one element.

The given eltype<:Number is the type of the returned vector.

# Example

```julia-repl
julia> nep=nep_gallery("dep0");
julia> x=ones(size(nep,1));
julia> s=compute_rf(Float64,nep,x)[1]; # Take just first element
0.6812131933795569
julia> x'*compute_Mlincomb(nep,s,x)
-8.881784197001252e-16
```
"""
    function compute_rf_new(T0::Type{T}, nep::NEP, x, inner_solver::InnerSolver;
                        y=x, target=zero(T), λ0=target,
                        TOL=eps(real(T))*1e3, max_iter=30,kwargs...) where T


        pnep=create_proj_NEP(nep);
        n=size(nep,1);
        set_projectmatrices!(pnep,reshape(y,n,1),reshape(x,n,1));

        local λv
        try
            (λv,xv)=inner_solve(inner_solver,T0,pnep,σ=target,λv=[λ0]);
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



function compute_rf_new(::Type{T},nep::NEP,x, inner_solver::PolyeigInnerSolver; y=x, target=zero(T), λ0=target,
                    TOL=eps(real(T))*1e3,max_iter=10) where T

    @show λ0
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
