module NEPSolver
    using ..NEPCore
    using ..NEPTypes
    using ..LinSolvers
    using LinearAlgebra
    using SparseArrays
    using Random

    export @ifd

    """
    @ifd(z)
Executes z if displaylevel>0.
"""
    macro ifd(z)
        return esc(:( if (displaylevel > 0); $z; end ))
    end

    export @ifdd

    """
    @ifd(z)
Executes z if displaylevel>1.
"""
    macro ifdd(z)
        return esc(:( if (displaylevel > 1); $z; end ))
    end

    ## NEP-Methods

    include("method_newton.jl")
    include("method_iar.jl")
    include("method_iar_chebyshev.jl")
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
    include("method_ilan.jl")
    include("method_ilan_benchmark.jl")
    include("method_nleigs.jl")

    include("inner_solver.jl");


    """
        σ = closest_to(λ_vec::Array{T,1},  λ::T) where {T<:Number}

    Finds the value `σ` in the vector `λ_vec` that is closes to the value `λ`.
    """
    function closest_to(λ_vec::Array{T,1}, λ::T) where {T<:Number}
        idx = argmin(abs.(λ_vec .- λ))
        return λ_vec[idx]
    end


end #End module
