export effenberger_deflation
export DeflatedNEP
export DeflatedNEPMM
export DeflatedGenericNEP
export DeflatedSPMF
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
     return DeflatedNEPMM(nep,S0,V0)
end


struct DeflatedNEPMM <: NEP
    orgnep::NEP
    S0
    V0
end

struct DeflatedSPMF{T,NEP1,NEP2} <: AbstractSPMF{T}
    orgnep::AbstractSPMF{T}
    spmf::SPMFSumNEP{NEP1,NEP2}
    S0
    V0
end

struct DeflatedGenericNEP <: NEP
    orgnep::NEP
    S0
    V0
end


DeflatedNEP=Union{DeflatedNEPMM,DeflatedGenericNEP,DeflatedSPMF{T} where T<:AbstractMatrix};

function size(nep::DeflatedNEP,dim=-1)
    n=size(nep.orgnep,1)+size(nep.V0,2);
    if (dim==-1)
        return (n,n)
    else
        return n
    end
end

## DeflatedGenericNEP optimized routines

function compute_Mder(nep::DeflatedGenericNEP,λ::Number,i::Integer=0)
end
function compute_Mlincomb(nep::DeflatedGenericNEP,λ::Number,
                 V::AbstractVecOrMat,a::Vector=ones(eltype(V),size(V,2)))
    # Do a factorization if c*(p*p*p+k*k*p*p) < k*k*p*p*p
    X=nep.V0;
    S=nep.S0;
    T=promote_type(typeof(λ),eltype(X),eltype(S));
    F=factorize(λ*I-S);
    k=size(V,2);
    Xhat=X/F;
    n0=size(nep.orgnep,1);

    p=size(S,1);
    # precompute terms with inverses:
    Q=Vector{Matrix}(undef,k)
    for i=0:k-1;
        QQ=zeros(T,p,k)
        QQ[:,i+1]=V[(n0+1):end,i+1];
        for j=(i-1):-1:0;
          QQ[:,j+1]=F\QQ[:,j+2];
          #QQ[:,j+1]=(λ*I-S)^(-(i-j))*V[(n0+1):end,i+1];
        end
        Q[i+1]= QQ
    end


    Z=zeros(eltype(V),n0,k);
    for j=0:k-1
        z=zeros(eltype(V),n0);
        for i=j:k-1
            #Vnew= (λ*I-S)^(-(i-j))*V[(n0+1):end,i+1]
            Vnew = Q[i+1][:,j+1];
            factor=((-1)^(i-j))*(a[i+1]*factorial(i)/factorial(j));
            z+=factor*Xhat*Vnew;
        end
        Z[:,j+1]=z;
    end
    Vnew=V[1:n0,:]*Diagonal(a)+Z;
    z_top=compute_Mlincomb(nep.orgnep,λ,Vnew);

    z_bottom=X'*V[1:n0,1]*a[1];
    return [z_top ; z_bottom];
end



## The DeflatedSPMF just delegates to the nep.spmf
compute_Mlincomb(nep::DeflatedSPMF,λ::Number,V::AbstractVecOrMat,a::Vector=ones(eltype(V),size(V,2)))= compute_Mlincomb(nep.spmf,λ,V,a);
compute_Mder(nep::DeflatedSPMF,λ::Number)=compute_Mder(nep.spmf,λ,0)
compute_Mder(nep::DeflatedSPMF,λ::Number,der)=compute_Mder(nep.spmf,λ,der)

## The DeflatedNEPMM calls corresponding MM-functions
function compute_MM(nep::DeflatedNEPMM,S,V)
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
# Use the MM to compute Mlincomb for DeflatedNEPMM
compute_Mlincomb(nep::DeflatedNEPMM,λ::Number,
                 V::AbstractVecOrMat,a::Vector=ones(eltype(V),size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)

function compute_Mder(nep::DeflatedNEPMM,λ::Number,i::Integer=0)
    # Use full to make it work with MSLP. This will not work for large and sparse.
    return Matrix(compute_Mder_from_MM(nep,λ,i));
end
