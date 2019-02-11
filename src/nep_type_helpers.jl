# Helper functions for NEPTypes (imported from NEPTypes.jl)
export Mder_NEP;
export Mder_Mlincomb_NEP;

# Type representing function handle NEP
struct Mder_Mlincomb_NEP <: NEP
    n::Int;
    Mder_fun::Function #
    maxder_Mder::Int;
    Mlincomb_fun::Function
    maxder_Mlincomb::Int;
end
"""
    Mder_NEP(n,Mder_fun; maxder=max)

Creates a `NEP` from its `compute_Mder` function defined by the
function handle `Mder_fun`. The `Mder_fun(λ,der)` takes two parameters a
scalar `λ::Number`, derivative number  `der`. The size `n::Int` must also
be specified. The function `Mder_fun(λ,der)` should return the n x n matrix
corresponding to the  `der`th derivatve. If only a
limited number of derivatives are available, `maxder` should be set, e.g.,
if not derivatives are implemented, set `maxder=0`. The function
`compute_Mlicomb` is automatically available by (delegation) to
`compute_Mlincomb_from_Mder`.

Note: This is a convenience function it is not recommended for high performance computations, since, e.g., it will not maintain type stability.

# Example
 The following example defines a linear eigenvalue problem
`A0+λA1` defined in an external function.
```julia-repl
julia> using LinearAlgebra; # For the I function
julia> function my_Mder(s,der)
    A0=ones(Float64,3,3)-I; A0[1,1]=-1;
    A1=ones(Float64,3,3)*3; A1[2,3]=0;
    if (der==0)
       return A0+A1*s;
    elseif (der==1)
       return A1;
    else
       return zero(A0);
    end
end
julia> nep=Mder_NEP(3,my_Mder);
julia> (λ,v)=augnewton(nep,v=ones(3));
julia> norm(compute_Mder(nep,λ)*v)
5.551115123125783e-17
```
"""
function Mder_NEP(n,Mder_fun::Function; maxder=typemax(Int64))
    maxder_Mlincomb=-1; # This means we delegate every Mlincomb_call
    Mlincomb_fun=s->s; # Dummy function
    return Mder_Mlincomb_NEP(n,
                             Mder_fun,maxder,
                             Mlincomb_fun,maxder_Mlincomb);
end


"""
    Mder_Mlincomb_NEP(n,Mder_fun, [maxder_Mder,] Mlincomb_fun, [maxder_Mlincomb])

Creates a `NEP` from its `compute_Mder` and `compute_Mlincomb`functions
 defined by the function handles `Mder_fun` and `Mlincomb_fun`. The `Mlincomb_fun(λ,X)` takes two parameters a
scalar `λ::Number` and a matrix `X`.  The size `n::Int` must also
be specified. The function `Mlincomb_fun(λ,X)` should return a vector corresponding of the linear combination of derivatives multiplied with the vectors in X. If only a
limited number of derivatives are implemented, `maxder_Mder` or
`maxder_Mlincomb` should be set, e.g., if not derivatives are implemented,
set `maxder=0`.

See also `Mder_NEP`.

Note: This is a convenience function it is not recommended for high performance computations, since, e.g., it will not maintain type stability.

# Example

The following example defines a linear eigenvalue problem
`A0+λA1` defined in an external function.
```julia-repl
julia> using LinearAlgebra; # For the I function
julia> function my_Mder(s,der)
    A0=ones(Float64,3,3)-I; A0[1,1]=-1;
    A1=ones(Float64,3,3)*3; A1[2,3]=0;
    if (der==0)
       return A0+A1*s;
    elseif (der==1)
       return A1;
    else
       return zero(A0);
    end
end
julia> function my_Mlincomb(s,X) # Compute linear comb of derivatives
    A0=ones(Float64,3,3)-I; A0[1,1]=-1;
    A1=ones(Float64,3,3)*3; A1[2,3]=0;
    if (size(X,2) <= 1)
       return A0*X[:,1]+s*A1*X[:,1]
    else # This means: size(X,2) => 2
       return A0*X[:,1]+A1*(s*X[:,1]+X[:,2]);
    end
end
julia> nep=Mder_Mlincomb_NEP(3,my_Mder,my_Mlincomb);
julia> (λ,v)=augnewton(nep,v=[1.0; 2.3; 0.0])
julia> norm(compute_Mder(nep,λ)*v)
6.798699777552591e-17
```
"""
function Mder_Mlincomb_NEP(n,
                           Mder_fun::Function, maxder_Mder,
                           Mlincomb_fun::Function, maxder_Mlincomb::Int=typemax(Int64))
    return Mder_Mlincomb_NEP(n,Mder_fun,maxder,
                             Mlincomb_fun,maxder_Mlincomb);
end
Mder_Mlincomb_NEP(n,Mder_fun::Function,Mlincomb_fun::Function) = Mder_Mlincomb_NEP(n,Mder_fun,typemax(Int64),Mlincomb_fun);



function compute_Mder(nep::Mder_Mlincomb_NEP,λ::Number,der::Integer=0)
    if (der>nep.maxder_Mder)
        error("Derivatives higher than ",nep.maxder_Mder," are not available.");
    end
    return nep.Mder_fun(λ,der);
end

function compute_Mlincomb(nep::Mder_Mlincomb_NEP,λ::Number,V::AbstractVecOrMat)
    if ((size(V,2)-1)<=nep.maxder_Mlincomb)
        return nep.Mlincomb_fun(λ,V); # Call external routine
    else
        # Delegate to if we do not have implementation of Mlincomb
        return compute_Mlincomb_from_Mder(nep,λ,V,ones(eltype(V),size(V,2)))
    end
end

function size(nep::Mder_Mlincomb_NEP)
    return (nep.n,nep.n)
end
function size(nep::Mder_Mlincomb_NEP,dim)
    return nep.n
end
