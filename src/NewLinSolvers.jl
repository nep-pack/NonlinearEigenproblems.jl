# imported from LinSolvers.jl

export LinSolverCreator, BackslashLinSolverCreator;
export FactorizeLinSolverCreator, DefaultLinSolverCreator;
export GMRESLinSolverCreator;

export create_linsolver;

abstract type LinSolverCreator ; end

struct BackslashLinSolverCreator <: LinSolverCreator
end

function create_linsolver(creator::BackslashLinSolverCreator,nep,λ)
    return BackslashLinSolver(nep,λ);
end


struct FactorizeLinSolverCreator <: LinSolverCreator
    umfpack_refinements::Int
    function FactorizeLinSolverCreator(umfpack_refinements::Int=1)
       return new(umfpack_refinements)
    end
end
# For the moment, Factorize is the default behaviour
DefaultLinSolverCreator = FactorizeLinSolverCreator


function create_linsolver(creator::FactorizeLinSolverCreator,nep,λ)
    return FactorizeLinSolver(nep,λ,creator.umfpack_refinements);
end

struct GMRESLinSolverCreator{T} <: LinSolverCreator where {T}
    kwargs::T
end
function GMRESLinSolverCreator(;kwargs...)
    return GMRESLinSolverCreator{typeof(kwargs)}(kwargs)
end


function create_linsolver(creator::GMRESLinSolverCreator,nep::NEP, λ)
    return GMRESLinSolver{typeof(λ)}(nep, λ, creator.kwargs)
end
