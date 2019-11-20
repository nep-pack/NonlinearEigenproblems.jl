module NEPTransformations

    # Transformations between types
    using SparseArrays

    include("nleigs_coefficients.jl")


    export taylor_expansion_pep
    export shift_and_scale
    export mobius_transform
    export CORKPencil
    export buildPencil
    export CorkLinearization
    export DefaultCorkLinearization
    export IarCorkLinearization
    export NleigsCorkLinearization
    export CORKPencilLR
    export lowRankCompress


    # We overload these
    import ..NEPCore.compute_Mder
    import ..NEPCore.compute_Mlincomb
    import ..NEPCore.compute_Mlincomb!
    import ..NEPCore.compute_MM
    import ..NEPCore.compute_resnorm
    import Base.size
    import SparseArrays.issparse

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

    # Example
    ```julia-repl
    julia> nep0=nep_gallery("pep0")
    julia> σ=3; α=10;
    julia> nep1=shift_and_scale(nep0,shift=σ,scale=α)
    julia> opnorm(compute_Mder(nep0,α*(4+4im)+σ)-compute_Mder(nep1,4+4im))
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
            AA = zero(At[j+1])
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
        J = Matrix{eltype(orgnep.A[1])}(I, size(orgnep.A[1]))
        return DEP([broadcast(*,orgnep.A,exp.(-orgnep.tauv*shift)/scale); [-shift/scale * J] ],[orgnep.tauv*scale; 0])
    end


    # Native implementation for SPMF. Only for generic SPMF (not AbstractSPMF)
    function shift_and_scale(orgnep::SPMF_NEP;shift=0,scale=1)
        orgfv=get_fv(orgnep);
        m=size(orgfv,1);
        fv = Vector{Function}(undef, m)
        # create anonymous functions corresponding to the
        # shift and scaled problem
        for i=1:m
            fv[i] = S -> orgfv[i](scale*S + shift*one(S))::((S isa Number) ? Number : Matrix)
        end
        # create a new SPMF with transformed functions
        return SPMF_NEP(get_Av(orgnep),fv;Schur_fact=orgnep.Schur_factorize_before);
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


    function compute_Mlincomb(nep::ShiftScaledNEP,λ::Number,V::AbstractVecOrMat,
                              a::Vector=ones(size(V,2)))
        # Multiply the coefficient matrix V with a diagonal matrix
        # corresponding to the chain rule
        p=size(V,2);
        z=nep.scale.^Array{eltype(V),1}(0:p-1)
        W=V * diagm(0 => z) # not very fast
        return compute_Mlincomb(nep.orgnep,nep.scale*λ+nep.shift,W,a);
    end



    """
        mobius_transform(orgnep::NEP,[,a=1][,b=0][,c=0][,d=1])

    Transforms a nep (orgnep) ``M(λ)v`` to a new nep ``T(λ)=M((aλ+b)/(cλ+d))``.
    This function tries to preserve the type such that `T`
    and `M` are of the same NEP-type (see `shift_and_scale()`).
    If it cannot be preserved it will return a `MobiusTransformedNEP`.
    It is in general advised to try to preserve the type, and
    the use of `MobiusTransformedNEP` can
    considerably slow down NEP-access.

    # Example

    ```julia-repl
    julia> nep0=nep_gallery("pep0")
    julia> a=1; b=3; c=4; d=5;
    julia> nep1=mobius_transform(nep0,a=a,b=b,c=c,d=d);
    julia> s=3;
    julia> opnorm(compute_Mder(nep0,(a*s+b)/(c*s+d))-compute_Mder(nep1,s))
    0.0
    ```
    """
    function mobius_transform(orgnep::NEP;a=1,b=0,c=0,d=1)
        return MobiusTransformedNEP(orgnep;a=a,b=b,c=c,d=d)
    end

    function mobius_transform(orgnep::SPMF_NEP;a=1,b=0,c=0,d=1)
        orgfv=get_fv(orgnep);
        m=size(orgfv,1);
        fv = Vector{Function}(undef, m)
        # create anonymous functions corresponding to the
        # möbius transformed problem
        for i=1:m
            fv[i] = S -> orgfv[i]((a*S + b*one(S)) / (c*S + d*one(S)))::((S isa Number) ? Number : Matrix)
        end
        # create a new SPMF with transformed functions
        return SPMF_NEP(get_Av(orgnep),fv;Schur_fact=orgnep.Schur_factorize_before);
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
        return compute_MM(nep.orgnep, (nep.a*S + nep.b*I) / (nep.c*S + nep.d*I), V)
    end

    function compute_Mlincomb(nep::MobiusTransformedNEP,λ::Number,V::AbstractVecOrMat
                              ,a::Vector=ones(size(V,2)))
        # I did not found a better way
        return compute_Mlincomb_from_MM(nep,λ,V,a)
    end

    # consider renaming this function
    """
        taylor_expansion_pep(orgnep::NEP[,d=2])

    Compute the truncated (with `d` terms) Taylor series of the NEP.
    The output is a [`PEP`](@ref).
    """
    function taylor_expansion_pep(nep::NEP,d::Integer=2)
        A=Vector{Matrix}(undef,d+1)
        for i=0:d
            A[i+1]=compute_Mder(nep,0,i)/factorial(i);
        end
        return PEP(A)
    end

"""
    struct CORKPencil
    function CORKPencil(M,N,Av,Bv,Z);
    function CORKPencil(nep,is)

The struct `CORKPencil` represents a pencil with a particular structure,
as given in the reference. The data can either be constructed directly
via the first constructor, or from a NEP in the second constructor.
The second constructor takes a `NEP` and  `is` which specifies
a CORK-structure as well as an approximation method. This can be objects
of the type `IarCorkLinearization` or `NleigsCorkLinearization`,
which are the CORK-linearizations equivalent to (certain versions of)
[`iar`](@ref) and [`nleigs`](@ref).

See [`buildPencil`](@ref) how to build standard pencil.


# Example:

The following example constructs a  `CORKPencil` from
a general NEP and then computes approximations of NEPs
by the interpolation approach of `nleigs`.

```julia-repl
julia> using LinearAlgebra
julia> A=(1:4)*(1:4)'+I; B=diagm(1 => [1,2,3]); C=ones(4,4);
julia> f1= λ-> one(λ);
julia> f2= λ-> λ;
julia> f3= λ-> exp(sin(λ/2));
julia> nep=SPMF_NEP([A,B,C],[f1,f2,f3]);
julia> cp=CORKPencil(nep,NleigsCorkLinearization());
julia> (A,B)=buildPencil(cp) # Form the pencil
julia> λv=eigen(A,B).values;
julia> λ=λv[sortperm(abs.(λv))[1]]; # Smallest eigval
julia> minimum(svdvals(compute_Mder(nep,λ))) # It's a solution
2.4364382475487156e-11
```

# References:

*  Compact rational Krylov methods for nonlinear eigenvalue problems SIAM Journal on Matrix Analysis and Applications, 36 (2), 820-838, 2015.

    """
struct CORKPencil{T1<:AbstractMatrix,T2<:AbstractMatrix}
        M::T1
        N::T1
        Av::Vector{T2}   # Array of Array of matrices
        Bv::Vector{T2}   # Array of Array of matrices
    end

    abstract type CorkLinearization end;

    struct DefaultCorkLinearization <: CorkLinearization end;
    # TODO: define default

    struct IarCorkLinearization <: CorkLinearization
        d::Int
        function IarCorkLinearization(;d=10)
            return new(d);
        end
    end

    struct NleigsCorkLinearization <: CorkLinearization
        Σ::Vector
        Ξ::Vector
        maxdgr::Int
        tollin::Float64;
        function NleigsCorkLinearization(;Σ=[-1.0-1im,-1+1im,+1+1im,1-1im],Ξ = [Inf], maxdgr=100 , tollin=1e-6 )
            return new(Σ,Ξ,maxdgr,tollin);
        end
    end

    function CORKPencil(nep,is::IarCorkLinearization)
        M=diagm( 0 =>  ones(is.d) )[2:end,:]
        N=diagm( -1 =>  1 ./ (1:is.d-1) )[2:end,:]
        Av=Array{AbstractMatrix,1}(undef, is.d)
        Bv=Array{AbstractMatrix,1}(undef, is.d)
        Av[1]=-compute_Mder(nep,0,0)
        for j=2:is.d Av[j]=zero(Av[1])              end
        for j=1:is.d Bv[j]=compute_Mder(nep,0,j)/j  end
        return CORKPencil(M,N,Av,Bv)
    end

    function CORKPencil(nep,is::NleigsCorkLinearization)
        D,β,ξ,σ=nleigs_coefficients(nep,is.Σ,tollin=is.tollin,Ξ=is.Ξ)
        d=length(β)-1
        σ=σ[1:d+1]; β=β[1:d+1]; ξ=ξ[1:d+1]
        σ,β,ξ=promote(σ,β,ξ) # be sure that M and N have the same type
        M=diagm( -1 => σ[1:d], 0 =>  β[1:d] )[2:end-1,1:end-1]
        N=diagm( -1 => ones(d), 0 =>  β[1:d]./ξ[1:d] )[2:end-1,1:end-1]
        Av=Array{AbstractMatrix,1}(undef, d)
        Av[1:d-1]=D[1:d-1]; Av[d]=D[d]-σ[d]/β[d+1]*D[d+1]
        Bv=Array{AbstractMatrix,1}(undef, d)
        Bv[1:d-1]=D[1:d-1]/ξ[d+1]; Bv[d]=D[d]/ξ[d+1]-D[d+1]/β[d+1]
        return CORKPencil(M,N,Av,Bv)
    end

    # construct the linearization
"""
    (A,B)=buildPencil(cp)

Constructs a pencil from a `CORKPencil` or `CORKPencilLR`. The returned
matrices correspond to a generalized eigenvalue problem. See also [`CORKPencil`](@ref), [`CORKPencilLR`](@ref).

# Example


"""
    function buildPencil(cp::CORKPencil)
        n=size(cp.Av[1],1)
        if issparse(cp.Av[1])   II=sparse(I, n, n)
        else                    II=Matrix(I, n, n)    end
        return [hcat(cp.Av...); kron(cp.M,II)], [hcat(cp.Bv...); kron(cp.N,II)]
    end

    """
        struct CORKPencilLR
        function CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z);

Represents / constructs a low-rank CORK-pencil. `AvLR`, `BvLR`
and `Z` correspond to the low-rank factorization of terms
`Av` and `Bv`. See [`CORKPencil`](@ref) and reference below.

    # Example:

The example illustrate a low-rank linearization of a Taylor expansion
of the NEP ``A0-λI+vv^Te^{-λ}``.

    ```julia-repl
    julia> A0=[1.0 3.0; -1.0 2.0]/10;
    julia> v=reshape([-1.0 ; 1]/sqrt(2),n,1);

    julia> Av=[-A0-v*v']
    julia> Bv=[-one(A0)-v*v']
    julia> BvLR=[v/2, -v/3, v/4, -v/5, v/6, -v/7,  v/8, -v/9]
    julia> AvLR=zero.(BvLR);
    julia> Z=v;
    julia> d=9;
    julia> M=diagm( 0 =>  ones(d) )[2:end,:]
    julia> N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]
    julia> cplr=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z);
    julia> (AA,BB)=build_CORKPencil(cplr2);
    julia> λ=eigen(AA,BB).values[end];
    julia> minimum(svdvals(A0-λ*I+v*v'*exp(-λ)))
    8.870165379112754e-13
    ```

# References:

*  Section 7 in Compact rational Krylov methods for nonlinear eigenvalue problems
SIAM Journal on Matrix Analysis and Applications, 36 (2), 820-838, 2015.
"""
    struct CORKPencilLR
    	M::AbstractMatrix
    	N::AbstractMatrix
    	Av::Vector{AbstractMatrix}
    	AvLR::Vector{AbstractMatrix} # In CORK-paper Avtilde
    	Bv::Vector{AbstractMatrix}
    	BvLR::Vector{AbstractMatrix} # In CORK-paper Bvtilde
        Z::AbstractMatrix
    end


"""
    cplr=lowRankCompress(cp_org::CORKPencil,dtilde,rk)

Constructs a `CORKPencilLR` from a `CORKPencil`. This is done by
assuming that terms higher than `dtilde` are of low rank, with rank `rk`.
More precisely, all `A[j]` and `B[j]` for `j>dtilde` are assumed
to be of the form ``C_jZ^T``.

# Example:

** TODO

"""
function lowRankCompress(cp_org::CORKPencil,dtilde,rk)
    d=length(cp_org.Av);
    # Only use the first low-rank term as a basis (this
    # can cause if the first low-rank term has lower rank
    # than the other ones).
    Z=svd(cp_org.Bv[dtilde+1]).V[:,1:rk];

    # Check that the M and N-matrices have the block triangular structure
    if ((norm(cp_org.M[1:(dtilde-1),(dtilde+1):end]) > 0) ||
        (norm(cp_org.N[1:(dtilde-1),(dtilde+1):end]) > 0))
        error("The M-matrix does not have the required structure. Try increasing dtilde.");
    end


    # Manually extract the low-rank coefficients by
    # multiplying by Z. Works since Z is orthogonal.
    Bvtilde=map(i-> cp_org.Bv[i]*Z,  (dtilde+1):d);
    Avtilde=map(i-> cp_org.Av[i]*Z,  (dtilde+1):d);

    return CORKPencilLR(cp_org.M,cp_org.N,
                        cp_org.Av[1:dtilde],Avtilde,
                        cp_org.Bv[1:dtilde],Bvtilde,
                        Z);
end
function buildPencil(cp::CORKPencilLR)

    # Extract the variables
    n=size(cp.Av[1],1)
    dtilde=length(cp.Av); # Nof full rank terms
    d=length(cp.Av)+length(cp.AvLR);
    rk=size(cp.Z,2);

    # Type logic. The pencil will be of the largest type of the eltypes
    T=promote_type(eltype(cp.Z),
                   eltype(cp.Av[1]),
                   eltype(cp.AvLR[1]),
                   eltype(cp.Bv[1]),
                   eltype(cp.BvLR[1]),
                   eltype(cp.M),
                   eltype(cp.N));


    In=Matrix{T}(I,n,n);
    Idtilde=Matrix{T}(I,rk,rk);

    ## Note: this info is not in CORK paper:
    ##    M11: (dtilde-1) x (dtilde)
    ##    M21: (d-dtilde) x (dtilde)
    ##    M22: (d-dtilde) x (d-dtilde)
    ## same for N

    # Extract Mij and Nij:
    M11=cp.M[1:(dtilde-1),1:dtilde];
    M21=cp.M[(dtilde):end,1:dtilde]
    M22=cp.M[(dtilde):end,(dtilde+1):end]

    N11=cp.N[1:(dtilde-1),1:dtilde];
    N21=cp.N[(dtilde):end,1:dtilde]
    N22=cp.N[(dtilde):end,(dtilde+1):end]

    # Build the B-matrix
    Btilde1=[hcat(cp.Bv[1:dtilde]...) hcat(cp.BvLR...) ]
    Btilde2=[kron(N11,In) zeros((dtilde-1)*n,(d-dtilde)*rk)]
    Btilde3=[kron(N21,cp.Z') kron(N22,Idtilde)];
    Btilde=vcat(Btilde1,Btilde2,Btilde3);


    # Build the A-matrix
    Atilde1=[hcat(cp.Av[1:dtilde]...) hcat(cp.AvLR...) ]
    Atilde2=[kron(M11,In) zeros((dtilde-1)*n,(d-dtilde)*rk)]
    Atilde3=[kron(M21,cp.Z') kron(M22,Idtilde)];
    Atilde=vcat(Atilde1,Atilde2,Atilde3);

    return (Atilde,Btilde);
end

end  # End Module
