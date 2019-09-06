module NEPTransformations

    # Transformations between types
    using SparseArrays

    include("nleigs_coefficients.jl")


    export transform_to_pep
    export shift_and_scale
    export mobius_transform
    export CORK_pencil
    export build_CORK_pencil
    export compute_CORK_pencil
    export CorkLinearization
    export DefaultCorkLinearization
    export IarCorkLinearization
    export NleigsCorkLinearization

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
        transform_to_pep(orgnep::NEP[,d=2])

    Compute the truncated (with `d` terms) Taylor series of the NEP.
    The output is a [`PEP`](@ref).
    """
    function transform_to_pep(nep::NEP,d::Integer=2)
        A=Array{Array{Float64, 2}}(d+1)
        for i=0:d
            A[i+1]=compute_Mder(nep,0,i)/factorial(i);
        end
        return PEP(A)
    end

    # representation of the structured linearizations used in the CORK framework
    struct CORK_pencil{T<:AbstractMatrix}
        M::T
        N::T
        Av::Vector{T}   # Array of Array of matrices
        Bv::Vector{T}   # Array of Array of matrices
    end

    abstract type CorkLinearization end;

    struct DefaultCorkLinearization <: CorkLinearization end;
    # define default

    struct IarCorkLinearization <: CorkLinearization
        d::Int
        function IarCorkLinearization(;d=10)
            return new(d);
        end
    end

    struct NleigsCorkLinearization <: CorkLinearization
        Σ::Union{Vector,Symbol}
        #Ξ::Union{Vector,Symbol} # add later
        maxdgr::Int
        tollin::Float64;
        function NleigsCorkLinearization(;Σ=[-1.0-1im,-1+1im,+1+1im,1-1im], maxdgr=100 , tollin=1e-6 )
            return new(Σ,maxdgr,tollin);
        end
    end

    function compute_CORK_pencil(nep,is::IarCorkLinearization)
        M=diagm( 0 =>  ones(is.d) )[2:end,:]
        N=diagm( -1 =>  1 ./ (1:is.d-1) )[2:end,:]
        Av=Array{AbstractMatrix,1}(undef, is.d)
        Bv=Array{AbstractMatrix,1}(undef, is.d)
        Av[1]=-compute_Mder(nep,0,0)
        for j=2:is.d Av[j]=zero(Av[1])              end
        for j=1:is.d Bv[j]=compute_Mder(nep,0,j)/j  end
        return CORK_pencil(M,N,Av,Bv)
    end

    function compute_CORK_pencil(nep,is::NleigsCorkLinearization)
        D,β,ξ,σ=nleigs_coefficients(nep,is.Σ,tollin=is.tollin)
        d=length(β)-1
        σ=σ[1:d+1]; β=β[1:d+1]; ξ=ξ[1:d+1]
        M=diagm( -1 => σ[1:d], 0 =>  β[1:d] )[2:end-1,1:end-1]
        N=diagm( -1 => ones(d), 0 =>  β[1:d]./ξ[1:d] )[2:end-1,1:end-1]
        Av=Array{AbstractMatrix,1}(undef, d)
        Av[1:d-1]=D[1:d-1]; Av[d]=D[d]-σ[d]/β[d+1]*D[d+1]
        Bv=Array{AbstractMatrix,1}(undef, d)
        Bv[1:d-1]=D[1:d-1]/ξ[d+1]; Bv[d]=D[d]/ξ[d+1]-D[d+1]/β[d+1]
        return CORK_pencil(M,N,Av,Bv)
    end

    # construct the linearization
    function build_CORK_pencil(cp::CORK_pencil)
        n=size(cp.Av[1],1)
        if issparse(cp.Av[1])   II=sparse(I, n, n)
        else                    II=Matrix(I, n, n)    end
        return [hcat(cp.Av...); kron(cp.M,II)], [hcat(cp.Bv...); kron(cp.N,II)]
    end

end  # End Module
