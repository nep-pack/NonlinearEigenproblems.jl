module NEPSolver
    using ..NEPCore
    using ..NEPTypes
    using ..LinSolvers

    export compute_eigvec_from_eigval_lu
    export compute_eigvec_from_eigval_lopcg

    export @ifd

    """
    @ifd(z)
Executes z if displaylevel>0.
"""
    macro ifd(z)
        return esc(:( if (displaylevel > 0); $z; end ))
    end

    ## NEP-Methods

    include("method_newton.jl")
    include("method_iar.jl")
    #include("method_iar_chebyshev.jl")
    include("method_tiar.jl")
    include("method_infbilanczos.jl")
    include("method_mslp.jl")
    include("method_companion.jl")
    include("method_nlar.jl")
    include("method_sgiter.jl")
    include("method_rfi.jl")
    include("method_jd.jl")
    include("method_beyncontour.jl")
    include("method_blocknewton.jl")
    include("method_broyden.jl")


    """
    ### compute_eigvec_from_eigval_lu
    Compute an eigenvector approximation from an accurate
    eigenvalue approximation with an LU approach. \\
    `nep` is the nonlinear eigenvalue problem \\
    `λ` is the accurate eigenvalue approximation \\
    `linsolvercreator` is the linsolver creator for M(λ) \\
    \\
    OBS: if a LinSolver `M0inv` for M(λ) exists, call this function as \\
    `compute_eigvec_from_eigval_lu(nep,λ,(nep, σ) -> M0inv)`
    """
    function compute_eigvec_from_eigval_lu(
        nep::NEP,
        λ::Number,
        linsolvercreator::Function
        )

        local M0inv::LinSolver = linsolvercreator(nep,λ);
        n=size(nep,1);

        if (isa(M0inv,DefaultLinSolver))
            F=M0inv.Afact
        else
            F=lufact(compute_Mder(nep,λ)); # This requires matrix access
        end
        x=[-F[:U][1:end-1,1:end-1]\F[:U][1:end-1,end]; 1]

        if issparse(nep)
            Q=zeros(Int64,n);   Q[F[:q]]=1:n;
            return x[Q]
        else
            return x
        end

    end


    """
    ### compute_eigvec_from_eigval_lopcg
    Compute an eigenvector approximation from an accurate
    eigenvalue approximation. \\
    This function uses the Locally optimal PCG (LOPCG) applied to the matrix
    M(λ) if the nep is symmetric, i.e., nep==nept, or to the matrix
    M(λ)^H M(λ). \\
    For a reference see
    Arbenz, Peter, Daniel Kressner, and D. M. E. Zürich.
    "Lecture notes on solving large scale eigenvalue problems."
    D-MATH, EHT Zurich 2 (2012). \\
    `nep` is the nonlinear eigenvalue problem \\
    `nept` is the transpose of the nonlinear eigenvalue problem \\
    `λ` is the accurate eigenvalue approximation \\
    """
    function compute_eigvec_from_eigval_lopcg(
        nep::NEP,
        nept::NEP,
        λ::Number;
        x=randn(size(nep,1)),
        tol=1e-6,
        maxit=10000,
        errmeasure::Function = default_errmeasure(nep::NEP)
        )

        A=v->compute_Mlincomb(nept,conj(λ),compute_Mlincomb(nep,λ,v))

        # initialization
        x/=norm(x); v=A(x); ρ=x⋅v; q=zeros(Complex128,size(nep,1));
        k=1; err=1; tol=1e-12
        while (k<maxit)&&(err>tol)
          g=v-ρ*x;
          aa=[x -g q]'*[v -A(g) A(q)]; aa=(aa+aa')/2;
          mm=[x -g q]'*[x -g q]; mm=(mm+mm')/2;
          D,V=eig(aa,mm); ii=indmin(abs.(D));
          ρ=D[ii]; δ=V[:,ii]; q=[-g q]*δ[2:end];
          x=δ[1]*x+q; x/=norm(x); v=A(x); k+=1
          err=errmeasure(λ,x)
        end
        return x
    end

    """
        σ = closest_to(λ_vec::Array{T,1},  λ::T) where {T<:Number}

    Finds the value `σ` in the vector `λ_vec` that is closes to the value `λ`.

    """
    function closest_to(λ_vec::Array{T,1},  λ::T) where {T<:Number}
        idx = indmin(abs.(λ_vec - λ))
        return λ_vec[idx]
    end


end #End module
