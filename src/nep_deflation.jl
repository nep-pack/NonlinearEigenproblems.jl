export effenberger_deflation
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
    if (nep isa AbstractSPMF)
        return DeflatedSPMF(nep,V0,S0)
    else
        return DeflatedGenericNEP(nep,V0,S0)
    end

end


struct DeflatedGenericNEP <: NEP
    orgnep::NEP
    V0
    S0
end

struct DeflatedSPMF{T} <: AbstractSPMF{T}
    orgnep::NEP
    V0::Matrix{T}
    S0::Matrix{T}
end


DeflatedNEP=Union{DeflatedGenericNEP,DeflatedSPMF{T} where T<:AbstractMatrix};

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
# Use the MM to compute Mlincomb for DeflatedNEP
compute_Mlincomb(nep::DeflatedNEP,λ::Number,
                 V::AbstractVecOrMat,a::Vector=ones(eltype(V),size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)

function compute_Mder(nep::DeflatedNEP,λ::Number,i::Integer=0)
    # Use full to make it work with MSLP. This will not work for large and sparse.
    return Matrix(compute_Mder_from_MM(nep,λ,i));
end
