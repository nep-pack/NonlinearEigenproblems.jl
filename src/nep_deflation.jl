export effenberger_deflation
export DeflatedNEP
export DeflatedNEPMM
export DeflatedGenericNEP
export DeflatedSPMF
export deflate_eigpair
export get_deflated_eigpairs


"""
    struct DeflatedNEPMM <: NEP

Represents a deflated NEP where the compute functions
are carried out via `compute_MM`. See further
documentation [`deflate_eigpair`](@ref).

"""
struct DeflatedNEPMM <: NEP
    orgnep::NEP
    S0
    V0
end


"""
    struct DeflatedSPMF <: AbstractSPMF

Represents a deflated NEP based on transforming
the deflated NEP to SPMF-form. See further
documentation [`deflate_eigpair`](@ref).
"""
struct DeflatedSPMF{T,NEP1,NEP2} <: AbstractSPMF{T}
    orgnep::AbstractSPMF{T}
    spmf::SPMFSumNEP{NEP1,NEP2}
    S0
    V0
end

"""
    struct DeflatedGenericNEP <: NEP

Represents a deflated NEP based on transforming
carrying out compute functions of derivatives by
binomial expansions. See further
documentation [`deflate_eigpair`](@ref).
"""
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

## DeflatedGenericNEP optimized routines based on binomial expansion for derivative
function compute_Mlincomb(nep::DeflatedGenericNEP,λ::Number,
                 V::AbstractVecOrMat,a::Vector=ones(eltype(V),size(V,2)))

    X=nep.V0;
    S=nep.S0;
    T=promote_type(typeof(λ),eltype(X),eltype(S),eltype(a),eltype(V));
    # It probably only makes sense to do a factorization
    # if c*(p*p*p+k*k*p*p) < k*k*p*p*p
    F=factorize(λ*I-S);
    k=size(V,2);
    Xhat=X/F;
    n0=size(nep.orgnep,1);

    p=size(S,1);
    # precompute terms with inverses:
    Q=Vector{Matrix{T}}(undef,k)
    for i=0:k-1;
        QQ=zeros(T,p,k)
        QQ[:,i+1]=V[(n0+1):end,i+1];
        for j=(i-1):-1:0;
          QQ[:,j+1]=F\QQ[:,j+2];
          #QQ[:,j+1]=(λ*I-S)^(-(i-j))*V[(n0+1):end,i+1];
        end
        Q[i+1]= QQ
    end

    Z=zeros(T,n0,k);
    for j=0:k-1
        z=zeros(T,n0);
        for i=j:k-1
            #Vnew= (λ*I-S)^(-(i-j))*V[(n0+1):end,i+1]w
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

compute_Mder(nep::DeflatedGenericNEP,λ::Number)=compute_Mder(nep,λ,0)
function compute_Mder(nep::DeflatedGenericNEP,λ::Number,der::Integer)
    X=nep.V0;
    S=nep.S0;
    T=promote_type(typeof(λ),eltype(X),eltype(S))

    n0=size(nep.orgnep,1);
    p=size(S,1);

    Q = deflated_nep_compute_Q(nep,λ,der)

    M0=compute_Mder(nep.orgnep,λ,der);
    if (M0 isa SparseMatrixCSC)
        # return a sparse matrix
        (II,JJ,VV)=findnz(M0);
        (II2,JJ2,VV2)=findnz(sparse(Q))
        JJ2=JJ2 .+  n0;
        II_merged=[II;II2];
        JJ_merged=[JJ;JJ2];
        VV_merged=[VV;VV2];

        # If zero'th derivative add (2,1)-block
        if (der==0)
            (II3,JJ3,VV3)=findnz(sparse(Matrix(X')));
            II_merged=[II_merged;II3 .+ n0];
            JJ_merged=[JJ_merged;JJ3];
            VV_merged=[VV_merged;VV3];
        end
        return sparse(II_merged,JJ_merged,VV_merged,n0+p,n0+p);
    else
        # Return a dense matrix
        if (der==0)
            return [M0 Q ; X'  zeros(T,p,p)];
        else
            return [M0 Q ; zeros(T,p,n0)  zeros(T,p,p)];
        end
    end
end


function deflated_nep_compute_Q(nep::DeflatedNEP,λ::Number,der::Integer)
    X=nep.V0
    S=nep.S0
    T=promote_type(typeof(λ),eltype(X),eltype(S))

    n0=size(nep.orgnep,1)
    p=size(S,1)

    Q=zeros(T,n0,p)
    F=factorize(λ*I-S)
    Vnew = X
    for i=der:-1:0
        # Equivalent to: Vnew= X*(λ*I-S)^(-(der-i+1))
        Vnew = Vnew/F
        factor=((-1)^(der-i))*(factorial(der)/factorial(i))
        # Equivalent to: Mi=compute_Mder(nep.orgnep,λ,i); Q+=Mi*(Vnew*factor);
        for j = 1:p
            Q[:,j] = Q[:,j] + compute_Mlincomb(nep.orgnep, λ, Vnew[:,j], [factor], i)
        end
    end
    return Q
end

## The DeflatedSPMF just delegates to the nep.spmf
compute_Mlincomb(nep::DeflatedSPMF,λ::Number,V::AbstractVecOrMat,a::Vector)= compute_Mlincomb(nep.spmf,λ,V,a);
compute_Mlincomb(nep::DeflatedSPMF,λ::Number,V::AbstractVecOrMat)= compute_Mlincomb(nep.spmf,λ,V);
compute_Mder(nep::DeflatedSPMF,λ::Number)=compute_Mder(nep.spmf,λ,0)
compute_Mder(nep::DeflatedSPMF,λ::Number,der::Integer)=compute_Mder(nep.spmf,λ,der)
compute_MM(nep::DeflatedSPMF,par...)=compute_MM(nep.spmf,par...)
get_Av(nep::DeflatedSPMF)=get_Av(nep.spmf)
get_fv(nep::DeflatedSPMF)=get_fv(nep.spmf)


## The DeflatedNEPMM & DeflatedGenericNEP calls corresponding MM-functions
function compute_MM(nep::Union{DeflatedNEPMM,DeflatedGenericNEP},S,V)
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
    return compute_Mder_from_MM(nep,λ,i);
end






# Creates and SPMF with deflated `S0` and `V0`.
function create_spmf_dnep(nep::AbstractSPMF,S0,V0)
    Av_org=get_Av(nep);
    fv_org=get_fv(nep);
    m=size(fv_org,1);
    p=size(V0,2);
    n0=size(nep,1);

    m1=m;     # size of "old" part
    m2=m*p+1; # size of "deflation" part

    # spmf1: Create the "old" part
    A1=Vector{eltype(Av_org)}(undef,m1);
    for k=1:m
        A0k=Av_org[k];
        if (eltype(A1) <:  SparseMatrixCSC)
            (II,JJ,VV)=findnz(A0k)
            A1[k]=sparse(II,JJ,VV,n0+p,n0+p);
        else
            A1[k]=zeros(eltype(A0k),n0+p,n0+p)
            A1[k][1:n0,1:n0]=A0k;
        end
    end
    spmf1=SPMF_NEP(A1,fv_org,check_consistency=false)
    # spmf2: Create the additional deflation terms:

    #    Promote rules for the eltype:
    #    We may need to increase the eltype type size, since S0,V0 can be complex
    T=promote_type(eltype(V0),eltype(S0),eltype(Av_org[1]));
    local T_LowRankFactor;
    if (eltype(Av_org) <: SparseMatrixCSC)
        T_LowRankFactor=SparseMatrixCSC{T,Int64};
    else
        T_LowRankFactor=Matrix{T};
    end
    L2=Vector{T_LowRankFactor}(undef,m2);
    U2=Vector{T_LowRankFactor}(undef,m2);
    fv2=Vector{Function}(undef,m2);
    (λtmp,X)=eigen(S0);
    λ::Vector{T}=λtmp[:]; # Ensure type
    count=0;
    for i=1:p
        y=(V0*(X[:,i]));
        ei=zeros(p); ei[i]=1;
        x=(ei'/X);
        for r=1:m
            count=count+1;
            # This will automatically convert to sparse / full
            L2[count] = reshape([(Av_org[r]*y) ;zeros(p)],n0+p,1);
            U2[count] = reshape([zeros(n0);x'],n0+p,1);
            fv2[count]=S-> (S-λ[i]*one(S))\fv_org[r](S);
        end
    end
    # The constant term
    L2[m*p+1]=[zeros(n0,p);Matrix{T}(I,p,p)]
    U2[m*p+1]=[Matrix(V0);zeros(p,p)]
    fv2[m*p+1]= S->one(S);
    spmf2=LowRankFactorizedNEP(L2,U2,fv2);

    return SumNEP(spmf1,spmf2);
end

"""
    normalize_schur_pair!(S,V)

Makes the Schur pair ``(S,V)`` in the sense that `V'*V` is the
identity matrix. This will not work if `size(V,2)>size(V,1)`.

"""
function normalize_schur_pair!(S,V)
    if (size(V,2)>size(V,1))
        @warn "Cannot normalize short and skinny V-matrices."
    else
        (QQ,RR)=qr(V);
        V[:]=Matrix(QQ); # Use skinny QR-factorization
        S[:]=(RR*S)/RR;
    end
end

# Determine and verify that the deflate mode is correct.
function verify_deflate_mode(nep::NEP,mode)
    if (mode==:Auto)
        if (nep isa AbstractSPMF)
            mode=:SPMF
        else
            mode=:Generic
        end
    end
    if ((mode == :SPMF) || (mode == :SPMFPlain)) &&
        !(nep isa AbstractSPMF)
        error("SPMF-mode only possible for `AbstractSPMF`-NEPs")
    end
    return mode;
end
function verify_deflate_mode(nep::DeflatedNEP,mode)
    if (nep isa DeflatedSPMF && ((mode == :SPMF) || (mode == :Auto)) )
       return :SPMF
    elseif (nep isa DeflatedNEPMM && ((mode == :MM) || (mode == :Auto)) )
       return :MM
    elseif (nep isa DeflatedGenericNEP && ((mode == :Generic) || (mode == :Auto)) )
       return :Generic
    else
       error("Unknown mode / type");
    end
end

"""
    dnep=deflate_eigpair(orgnep::NEP,λ,v,[mode=:Auto])

This function creates a deflated NEP based on ``(λ,v)``, which
are assumed to an eigenpair of `nep`. Effectively,
the function will return `dnep::NEP` which has the same
solutions as orgnep, except those corresponding to ``(λ,v)``.
Deflation is typically used to avoid reconvergence.

If `orgnep` is a `DeflatedNEP`, the `orgnep` the
deflation in `orgnep` will be updated.

The `mode` kwarg can be `:Auto`, `:Generic`, `:SPMF`,
`:MM`. This specifies how the deflated NEP should be represented.
Which mode is the most efficient depends on
many problem properties.
If the original NEP is an `AbstractSPMF` with only
a few terms, `mode=:SPMF` may be efficient. The SPMF-mode
is based on a diagonalization of the deflated invariant
pair and is not necessarily robust when you deflate eigenvalues
near to each other. When `mode=:MM`
is used, all compute functions are implemented via calls
to the `compute_MM`. This can work well for small dense
problems. The `:Generic` is based on an explicit derivation
of the problem (via binomial expansion) which can be
efficient if low order derivates are needed. If
`:Auto` is selected, NEP-PACK tries to determine which
one is the most efficient based on the `orgnep`.

# Example:
```julia-repl
julia> nep=nep_gallery("dep0");
julia> (λ,v)=newton(nep);
julia> n=size(nep,1);
julia> dnep=deflate_eigpair(nep,λ,v)
julia> (λ2,v2)=augnewton(dnep);  # this converges to different eigval
julia> using LinearAlgebra;
julia> minimum(svdvals(compute_Mder(nep,λ2)))
9.323003321058995e-17
```
The function [`get_deflated_eigpairs()`](@ref) extracts the eigenpairs that
have been deflated. The returned pairs are eigenpairs of
the original NEP:
```
julia> dnep=deflate_eigpair(dnep,λ2,v2)
julia> (D,V)=get_deflated_eigpairs(dnep)
julia> norm(compute_Mlincomb(nep,D[1],V[:,1]))
6.164690797405912e-16
julia> norm(compute_Mlincomb(nep,D[2],V[:,2]))
5.20740757162067e-16
```

# References
* C. Effenberger, Robust solution methods for nonlinear eigenvalue problems, PhD thesis, 2013, EPF Lausanne
"""
function deflate_eigpair(nep::NEP,λ,v;mode=:Auto)
    mode=verify_deflate_mode(nep,mode);
    n=size(nep,1);
    S0=reshape([λ],1,1);
    V0=reshape(v,n,1);
    normalize_schur_pair!(S0,V0);
    # Create it based on nep
    return internal_create_deflated_nep(nep,S0,V0,mode);
end

function deflate_eigpair(nep::DeflatedNEP,λ,v;mode=:Auto)
    mode=verify_deflate_mode(nep,mode)

    T=promote_type(typeof(λ),eltype(v),eltype(nep.V0),eltype(nep.S0));

    n=size(nep.orgnep,1);
    p0=size(nep.V0,2);
    # fetch pair + expand with λ v
    V1=zeros(T,n,p0+1);
    S1=zeros(T,p0+1,p0+1);
    V1[1:n,1:end-1]=nep.V0
    V1[1:n,end]=v[1:n];
    S1[1:end-1,1:end-1]=nep.S0;
    S1[1:end,end]=[v[n+1:end];λ];

    # normalize schur pair
    normalize_schur_pair!(S1,V1);
    # Create it based on nep.orgnep
    return internal_create_deflated_nep(nep.orgnep,S1,V1,mode);
end

# Internal helper function to create a nep based on mode
function internal_create_deflated_nep(nep,S1,V1,mode)
    # create new DeflatedNEP
    if (mode==:MM)
        newnep=DeflatedNEPMM(nep,S1,V1);
        return newnep;
    elseif (mode==:SPMF)
        spmf=create_spmf_dnep(nep,S1,V1);;
        newnep=DeflatedSPMF(nep,spmf,S1,V1);
        return newnep;
    elseif (mode==:Generic)
        newnep=DeflatedGenericNEP(nep,S1,V1);
        return newnep;
    end
end


"""
    (D,V)=get_deflated_eigpairs(dnep::DeflatedNEP [λ,v])

Returns a vector of eigenvalues `D` and a matrix with corresponding
eigenvectors `V`. The eigenpairs correspond to the
original problem, underlying the `DeflatedNEP`.
The optional parameters `λ,v` allows the inclusion
of an additional eigpair. Essentially, the optional parameters
are the expanding the deflation and the running `get_deflated_eigpairs`
 with kwarg, i.e.,
```julia
(D,V)=get_deflated_eigpairs(deflate_eigpair(dnep,λ,v))`
```
See example in [`deflate_eigpair`](@ref).
"""
function get_deflated_eigpairs(nep::DeflatedNEP)
   V=nep.V0;
   S=nep.S0;
   (D,X)=eigen(S);
   return D,V[1:size(nep.orgnep,1),:]*X;
end

function get_deflated_eigpairs(nep::DeflatedNEP,λ,v)
   T=promote_type(typeof(λ),eltype(v),eltype(nep.V0),eltype(nep.S0));
   p0=size(nep.V0,2);
   n=size(nep.orgnep,1);
   # fetch pair + expand with λ v
   V1=zeros(T,n,p0+1);
   S1=zeros(T,p0+1,p0+1);
   V1[1:n,1:end-1]=nep.V0
   V1[1:n,end]=v[1:n];
   S1[1:end-1,1:end-1]=nep.S0;
   S1[1:end,end]=[v[n+1:end];λ];
   V=V1
   S=S1;
   (D,X)=eigen(S);
   return D,V[1:size(nep.orgnep,1),:]*X;
end
