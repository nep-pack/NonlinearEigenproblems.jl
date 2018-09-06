using LinearAlgebra
using SuiteSparse

struct LUFactors{T}
    A::AbstractMatrix
    LU::T

    # only store A for full matrices (needed to improve the solution accuracy)
    LUFactors(A::U, LU::T) where {U<:AbstractMatrix, T<:LU} = new{T}(A, LU)
    LUFactors(A::U, LU::T) where {U<:AbstractMatrix, T<:SuiteSparse.UMFPACK.UmfpackLU} = new{T}(reshape([], (0, 0)), LU)
end

struct LUCache{T<:Number}
    lu::Dict{T,LUFactors}
    verbose::Bool
end

LUCache(::Type{T}, verbose::Bool) where T<:Number = LUCache(Dict{T,LUFactors}(), verbose)

function solve(LU::LUFactors{LU{Complex{Float64},Matrix{Complex{Float64}}}}, y::AbstractVector{<:Number})
    x = LU.LU \ y
    # improve accuracy
    resid = y - LU.A * x
    err = LU.LU \ resid
    x + err
end

function solve(LU::LUFactors{SuiteSparse.UMFPACK.UmfpackLU{Complex{Float64},Int64}}, y::AbstractVector{<:Number})
    LU.LU \ y
    # improving accuracy is not needed with UmfpackLU
end

function lufactors(fun, v, verbose)
    verbose && println("Cache miss, computing LU factors for $v")
    A = fun(v)
    LUFactors(A, lu(A))
end

"Solve x for fun(v)*x = y, with cached lu-factors."
function lusolve(cache, fun, v, y)
    cache.verbose && println("LU solve for $v (cache size: $(length(cache.lu)))")
    lu = get!(() -> lufactors(fun, v, cache.verbose), cache.lu, v)
    solve(lu, y)
end
