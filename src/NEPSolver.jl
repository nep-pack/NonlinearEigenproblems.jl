module NEPSolver
    using NEPCore
    using NEPTypes
    using LinSolvers

    export compute_eigvec_from_eigval
    export @ifd

    """
    @ifd(z)
Executes z if displaylevel>0.
"""
    macro ifd(z)
        return :( if (displaylevel>0); $z; end )
    end

    ## NEP-Methods

    include("method_newton.jl")
    include("method_iar.jl")
    include("method_tiar.jl")
    include("method_infbilanczos.jl")
    include("method_mslp.jl")
    include("method_companion.jl")
    include("method_nlar.jl")
    include("method_sg.jl")
    include("method_rfi.jl")
#    include("method_jd_lin.jl")
    include("method_jd_quad.jl")


    """
    ### compute_eigvec_from_eigval
    Compute an eigenvector approximation from an accurate
    eigenvalue approximation. \\
    `nep` is the nonlinear eigenvalue problem \\
    `λ` is the accurate eigenvalue approximation \\
    `linsolvercreator` is the linsolver creator for M(λ) \\
    \\
    OBS: if a LinSolver `M0inv` for M(λ) exists, call this function as \\
    `compute_eigvec_from_eigval(nep,λ,(nep, σ) -> M0inv)`
    """
    function compute_eigvec_from_eigval(nep::NEP,λ,linsolvercreator::Function)

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




end #End module
