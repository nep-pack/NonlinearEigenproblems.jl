# Transformations between types
export ShiftScaleNEP
export MobiusTransformNEP
export transform_to_pep



   """
    type ShiftScaleNEP <: NEP
    ShiftScaleNEP(orgnep::NEP[,shift=0][,scale=1])

Transforms a nep (orgnep) M(λ)v to a new nep T(λ)=M(scale*λ+shift). This can be used if the method does not have an easy implementation of shift and scaling. Usage of this transformation can slow down the algorithm.

"""

type ShiftScaleNEP <: NEP
    shift::Number
    scale::Number
    orgnep::NEP
    function ShiftScaleNEP(orgnep;shift=0,scale=1)
        return new(shift,scale,orgnep)
    end
end



# Size and issparse implemented by delegation
function size(nep::ShiftScaleNEP,dim=-1)
    return size(nep.orgnep,dim)
end
function issparse(nep::ShiftScaleNEP)
    return issparse(nep.orgnep)
end


function compute_Mder(nep::ShiftScaleNEP,λ::Number,i::Integer=0)
    # Just scale the orgnep call. Chain rule for differentiation.
    return (nep.scale^i)*compute_Mder(nep.orgnep,nep.scale*λ+nep.shift,i);
end

function compute_MM(nep::ShiftScaleNEP,S,V)
    # Just call orgnep with a different S
    return compute_MM(nep.orgnep,S*nep.scale+eye(S)*nep.shift,V)
end


function compute_Mlincomb(nep::ShiftScaleNEP,λ::Number,V;a=ones(size(V,2)))
    # Multiply the coefficient matrix V with a diagonal matrix
    # corresponding to the chain rule
    p=size(V,2);
    z=nep.scale.^Array{eltype(V),1}(0:p-1)
    W=V*diagm(z); # not very fast
    return compute_Mlincomb(nep.orgnep,nep.scale*λ+nep.shift,W,a=a);
end






   """
    type MobiusTransformNEP <: NEP
    MobiusTransformNEP(orgnep::NEP[,a=1][,b=0][,c=0][,d=1])

Transforms a nep (orgnep) M(λ)v to a new nep T(λ)=M((a*λ+b)/(c*λ+d)). Usage of this transformation can slow down the algorithm.

"""

type MobiusTransformNEP <: NEP
    a::Number
    b::Number
    c::Number
    d::Number
    orgnep::NEP
    function MobiusTransformNEP(orgnep;a=1,b=0,c=0,d=1)
        return new(a,b,c,d,orgnep)
    end
end

# Size and issparse implemented by delegation
function size(nep::MobiusTransformNEP,dim=-1)
    return size(nep.orgnep,dim)
end
function issparse(nep::MobiusTransformNEP)
    return issparse(nep.orgnep)
end

function compute_Mder(nep::MobiusTransformNEP,λ::Number,i::Integer=0)
    # I did not found a better way
    return compute_Mder_from_MM(nep,λ,i);
end

function compute_MM(nep::MobiusTransformNEP,S,V)
    # Just call orgnep with a different S
    return compute_MM(nep.orgnep,(nep.a*S+nep.b*eye(S))/(nep.c*S+nep.d*eye(S)),V)
end

function compute_Mlincomb(nep::MobiusTransformNEP,λ::Number,V;a=ones(size(V,2)))
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
