macro log(level, string)
#    :( if verbose > $level println($string) end )
    :( if 2*0 > $level println($(esc(string))) end ) # TODO log level
end

struct LUFactors{T}
    A::AbstractArray
    LU::T

    # only store A for full matrices (needed to improve the solution accuracy)
    LUFactors(A::U, LU::T) where {U<:AbstractArray, T<:Base.LinAlg.LU} = new{T}(A, LU)
#    LUFactors(A::U, LU::T) where {U<:AbstractArray, T<:Base.SparseArrays.UMFPACK.UmfpackLU} = new{T}(A, LU)
    LUFactors(A::U, LU::T) where {U<:AbstractArray, T<:Base.SparseArrays.UMFPACK.UmfpackLU} = new{T}([], LU)
end

function solve(LU::LUFactors{Base.LinAlg.LU{Complex{Float64},Matrix{Complex{Float64}}}}, y::AbstractVector{<:Number})
    x = LU.LU \ y
    # improve accuracy
    resid = y - LU.A * x
    err = LU.LU \ resid
    x + err
end

function solve(LU::LUFactors{Base.SparseArrays.UMFPACK.UmfpackLU{Complex{Float64},Int64}}, y::AbstractVector{<:Number})
    x = LU.LU \ y
    # improve accuracy; seems not needed with UmfpackLU
#    resid = y - LU.A * x
#    err = LU.LU \ resid
#    x + err
end

lucache = Dict{Number,LUFactors}()

function lureset()
    @log(1, "Reset LU factor cache (cache size: $(length(lucache)))")
    global lucache = Dict{Number,LUFactors}()
end

function lufactors(funA, v)
    @log(1, "Cache miss, computing LU factors for $v")
    A = funA(v)
    LUFactors(A, lufact(A))
end

# lusolve: Solve x for funA(v)*x = y, with cached lu-factors
#   v  value for evaluating funA(v)
#   y  right hand side
function lusolve(funA, v::T, y::Vector{T}) where T<:Number
    @log(1, "LU solve for $v (cache size: $(length(lucache)))")
    lu = get!(() -> lufactors(funA, v), lucache, v)
    solve(lu, y)
end
