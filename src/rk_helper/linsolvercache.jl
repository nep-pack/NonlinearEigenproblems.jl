using NonlinearEigenproblems
using LinearAlgebra
using SuiteSparse

export LinSolverCache, solve

struct LinSolverCache{T<:Number}
    solver::Dict{T,LinSolver}
    nep::NEP
    linsolvercreator::Function
end

LinSolverCache(::Type{T}, nep, linsolvercreator) where T<:Number =
    LinSolverCache(Dict{T,LinSolver}(), nep, linsolvercreator)

function solve(cache, σ, y, add_to_cache)
    @debug "Lin solve for $σ, cache=$add_to_cache, size=$(length(cache.solver))"
    if add_to_cache
        solver = get!(() -> (
            @debug "Cache miss, creating lin solver for $σ";
            cache.linsolvercreator(cache.nep, σ)), cache.solver, σ)
    else
        solver = cache.linsolvercreator(cache.nep, σ)
    end
    lin_solve(solver, y)
end
