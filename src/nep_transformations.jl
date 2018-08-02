
# Transformations between types
export transform_to_pep
export shift_and_scale
export mobius_transform
export effenberger_deflation


"""
    struct ShiftScaledNEP <: NEP
    ShiftScaleNEP(orgnep::NEP[,shift=0][,scale=1])

Transforms a nep (orgnep) M(λ)v to a new nep T(λ)=M(scale*λ+shift). This can be used if the method does not have an easy implementation of shift and scaling. Usage of this transformation can slow down the algorithm.

"""

struct ShiftScaledNEP <: NEP
    shift::Number
    scale::Number
    orgnep::NEP
    function ShiftScaledNEP(orgnep;shift=0,scale=1)
        return new(shift,scale,orgnep)
    end
end
"""
    shift_and_scale(orgnep::NEP;shift=0,scale=1)
Transforms the orgnep by defining a new NEP from the relation
T(λ)=M(scale * λ+shift) where M is the orgnep. This function tries
 to preserve the NEP type, e.g., a shift_and_scale operation on
an SPMF-object, return an SPMF object. If it cannot preserve
the type, it will return a nep of the struct `ShiftScaledNEP`.

#    Example
```julia-repl
julia> nep0=nep_gallery("pep0")
julia> σ=3; α=10;
julia> nep1=shift_and_scale(nep0,shift=σ,scale=α)
julia> norm(compute_Mder(nep0,α*(4+4im)+σ)-compute_Mder(nep1,4+4im))
8.875435870738592e-12
```
"""
function shift_and_scale(orgnep::NEP;shift=0,scale=1)
    return ShiftScaledNEP(orgnep,shift=shift,scale=scale);
end

# Native implementation for PEPs by transformation of the coefficients. This is not necessarily fast.
function shift_and_scale(orgnep::PEP;shift=0,scale=1)
    At=copy(shift*scale*orgnep.A); # Allocate correct type (not so fast)
    m=size(At,1)-1
    for j=0:m
        AA=zeros(At[j+1])
        for i=j:m
            factor=((scale^j)*(shift^(i-j))*factorial(i)/factorial(i-j))/factorial(j)
            AA+=orgnep.A[i+1]*factor;


        end
        At[j+1]=AA;
    end
    return PEP(At)
end

# A shift and rescale of a DEP is DEP where the matrix coefficients and the delays are rescaled and another matrix coefficient with delay zero is added.
function shift_and_scale(orgnep::DEP;shift=0,scale=1)
    return DEP([broadcast(*,orgnep.A,exp.(-orgnep.tauv*shift)/scale); [-shift/scale*eye(orgnep.A[1])] ],[orgnep.tauv*scale; 0])
end


# Native implementation for SPMF. Only for generic SPMF (not AbstractSPMF)
function shift_and_scale(orgnep::SPMF_NEP;shift=0,scale=1)
    orgfv=get_fv(orgnep);
    m=size(orgfv,1);
    fv=Array{Function}(m)
    # create anonymous functions corresponding to the
    # shift and scaled problem
    for i=1:m
        fv[i]= S -> orgfv[i](scale*S+shift*eye(S))
    end
    # create a new SPMF with transformed functions
    return SPMF_NEP(get_Av(orgnep),fv,orgnep.Schur_factorize_before);
end



# Size and issparse implemented by delegation
function size(nep::ShiftScaledNEP,dim=-1)
    return size(nep.orgnep,dim)
end
function issparse(nep::ShiftScaledNEP)
    return issparse(nep.orgnep)
end


function compute_Mder(nep::ShiftScaledNEP,λ::Number,i::Integer=0)
    # Just scale the orgnep call. Chain rule for differentiation.
    return (nep.scale^i)*compute_Mder(nep.orgnep,nep.scale*λ+nep.shift,i);
end

function compute_MM(nep::ShiftScaledNEP,S,V)
    # Just call orgnep with a different S
    return compute_MM(nep.orgnep,S*nep.scale+eye(S)*nep.shift,V)
end


function compute_Mlincomb(nep::ShiftScaledNEP,λ::Number,V;a=ones(size(V,2)))
    # Multiply the coefficient matrix V with a diagonal matrix
    # corresponding to the chain rule
    p=size(V,2);
    z=nep.scale.^Array{eltype(V),1}(0:p-1)
    W=V*diagm(z); # not very fast
    return compute_Mlincomb(nep.orgnep,nep.scale*λ+nep.shift,W,a=a);
end



"""
    mobius_transform(orgnep::NEP,[,a=1][,b=0][,c=0][,d=1])
Transforms a nep (orgnep) M(λ)v to a new nep T(λ)=M((a*λ+b)/(c*λ+d)).
This function tries to preserve the type such that T
and M are of the same NEP-type (see `shift_and_scale()`).
If it cannot be preserved it will return a `MobiusTransformedNEP`.
The use of `MobiusTransformedNEP` can considerably slow down the algorithm.

# Example
julia> nep0=nep_gallery("pep0")
julia> a=1; b=3; c=4; d=5;
julia> nep1=mobius_transform(nep0,a=a,b=b,c=c,d=d);
julia> s=3;
julia> norm(compute_Mder(nep0,(a*s+b)/(c*s+d))-compute_Mder(nep1,s))
0.0
"""

function mobius_transform(orgnep::NEP;a=1,b=0,c=0,d=1)
    return MobiusTransformedNEP(orgnep;a=a,b=b,c=c,d=d)
end

function mobius_transform(orgnep::SPMF_NEP;a=1,b=0,c=0,d=1)
    orgfv=get_fv(orgnep);
    m=size(orgfv,1);
    fv=Array{Function}(m)
    # create anonymous functions corresponding to the
    # möbius transformed problem
    for i=1:m
        fv[i]= S -> orgfv[i]((a*S+b*eye(S))/(c*S+d*eye(S)))
    end
    # create a new SPMF with transformed functions
    return SPMF_NEP(get_Av(orgnep),fv,orgnep.Schur_factorize_before);
end



"""
Hidden type representing a transformed NEP
"""

struct MobiusTransformedNEP <: NEP
    a::Number
    b::Number
    c::Number
    d::Number
    orgnep::NEP
    function MobiusTransformedNEP(orgnep;a=1,b=0,c=0,d=1)
        return new(a,b,c,d,orgnep)
    end
end

# Size and issparse implemented by delegation
function size(nep::MobiusTransformedNEP,dim=-1)
    return size(nep.orgnep,dim)
end
function issparse(nep::MobiusTransformedNEP)
    return issparse(nep.orgnep)
end

function compute_Mder(nep::MobiusTransformedNEP,λ::Number,i::Integer=0)
    # I did not found a better way
    return compute_Mder_from_MM(nep,λ,i);
end

function compute_MM(nep::MobiusTransformedNEP,S,V)
    # Just call orgnep with a different S
    return compute_MM(nep.orgnep,(nep.a*S+nep.b*eye(S))/(nep.c*S+nep.d*eye(S)),V)
end

function compute_Mlincomb(nep::MobiusTransformedNEP,λ::Number,V;a=ones(size(V,2)))
    # I did not found a better way
    return compute_Mlincomb_from_MM(nep,λ,V,a)
end


   """
   transform_to_pep(orgnep::NEP[,d=2])

Compute the truncated (with d term) Taylor series of a nep. The output is a PEP.
"""
# consider rename this function
function transform_to_pep(nep::NEP,d::Integer=2)
    A=Array{Array{Float64, 2}}(d+1)
    for i=0:d
        A[i+1]=compute_Mder(nep,0,i)/factorial(i);
    end
    return PEP(A)
end



"""
    effenberger_deflation(nep::NEP,S0,V0)

This function creates a deflated NEP based on (S0,V0), which
are assumed to an invariant pair of `nep`. Effectively,
the function should return a NEP which has the same
solutions as orgnep, except those corresponding to (S0,V0).


# Example:
```julia-repl
julia> nep=nep_gallery("dep0");
julia> (λ,v)=newton(nep);
julia> n=size(nep,1);
julia> S0=reshape([λ],1,1);
julia> V0=reshape(v,n,1);
julia> dnep=effenberger_deflation(nep,S0,V0)
julia> (λ2,v2)=augnewton(dnep);  # this converges to different eigval
julia> minimum(svdvals(compute_Mder(nep,λ2)))
9.323003321058995e-17
```

# References
* C. Effenberger, Robust solution methods for nonlinear eigenvalue problems, PhD thesis, 2013, EPF Lausanne 

"""

function effenberger_deflation(nep::NEP,S0,V0)
    return DeflatedNEP(nep,V0,S0)
end


struct DeflatedNEP <: NEP
    orgnep::NEP
    V0
    S0
end

function size(nep::DeflatedNEP,dim=-1)
    n=size(nep.orgnep,1)+size(nep.V0,2);
    if (dim==-1)
        return (n,n)
    else
        return n
    end
end

function compute_MM(nep::DeflatedNEP,S,V)
    orgnep=nep.orgnep;
    n0=size(orgnep,1);
    S0=nep.S0; V0=nep.V0
    p0=size(S0,1); p=size(S,1);
    V1=view(V,1:n0,            1:size(V,2))
    V2=view(V,(n0+1):size(V,1),1:size(V,2))
    Stilde=hcat(vcat(S0,zeros(p,p0)),vcat(V2,S))
    Vtilde=hcat(V0,V1);
    R=compute_MM(orgnep, Stilde,Vtilde);
    return vcat(R[1:n0,(size(nep.S0,1)+1):end],V0'*V1);
end

function compute_Mder(nep::DeflatedNEP,λ::Number,i::Integer=0)
    # Use full to make it work with MSLP. This will not work for large and sparse.
    return Matrix(compute_Mder_from_MM(nep,λ,i));
end
