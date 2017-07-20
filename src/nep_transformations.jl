# Transformations between types
export ShiftScaleNEP
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






