# imported from LinSolvers.jl

export LinSolverCreator, BackslashLinSolverCreator;
export FactorizeLinSolverCreator, DefaultLinSolverCreator;
export create_linsolver;
abstract type LinSolverCreator ; end

struct BackslashLinSolverCreator <: LinSolverCreator
end

function create_linsolver(creator::BackslashLinSolverCreator,nep,λ)
    return BackslashLinSolver(nep,λ);
end


struct FactorizeLinSolverCreator <: LinSolverCreator
    umfpack_refinements::Int
end
# For the moment, Factorize is the default behaviour
DefaultLinSolverCreator = FactorizeLinSolverCreator

function FactorizeLinSolverCreator(umfpack_refinements=1)
    return FactorizeLinSolverCreator(umfpack_refinements)
end


function create_linsolver(creator::FactorizeLinSolverCreator,nep,λ)
    return DefaultLinSolver(nep,λ,creator.umfpack_refinements);
end
