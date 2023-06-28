# imported from LinSolvers.jl

export LinSolverCreator, BackslashLinSolverCreator;
export FactorizeLinSolverCreator, DefaultLinSolverCreator;
export GMRESLinSolverCreator;
export DeflatedNEPLinSolverCreator

export create_linsolver;
using ..NEPTypes

abstract type LinSolverCreator ; end

"""
    struct BackslashLinSolverCreator <: LinSolverCreator

Creator to for the `BackslashLinSolver`, i.e., usage of backslash to make linear solves.
Specify objects of this type if you want the solver to use backslash.

See also: [`BackslashLinSolver`](@ref), [`create_linsolver`](@ref)
"""
struct BackslashLinSolverCreator <: LinSolverCreator
end

"""
    create_linsolver(creator::LinSovlerCreator,nep,λ)

Creates a `LinSolver` instance for the `nep` corresponding
which is evaluated in `λ`. The type of the output is
decided by dispatch and the type of the `LinSolverCreator`.

See also: `LinSolver`, `FactorizeLinSolverCreator`,
`BackslashLinSolverCreator`, `DefaultLinSolverCreator`,
`GMRESLinSolverCreator`.
"""
function create_linsolver(creator::BackslashLinSolverCreator,nep,λ)
    return BackslashLinSolver(nep,λ);
end

"""
    FactorizeLinSolverCreator(;umfpack_refinements,max_factorizations,nep,precomp_values)

`FactorizeLinSolverCreator`-objects can instantiate `FactorizeLinSolver`
objects via the `create_linsolver` function.

The `FactorizeLinSolver` is based on `factorize`-calls.
The time point of the call to `factorize` can be controlled by parameters
to `FactorizeLinSolverCreator`:

* By default, the `factorize` call is carried out by the instantiation of the `FactorizeLinSolver`, i.e., when the NEP-solver calls `create_linsolver`.

* You can also precompute the factorization, at the time point when you instantiate `FactorizeLinSolverCreator`. If you set `precomp_values::Vector{Number}` to a non-empty vector, and set `nep` kwarg, the factorization (of all λ-values in the `precomp_values`) will be computed  when the `FactorizeLinSolverCreator` is instantiated. If the NEP-solver calls a `create_linsolver` with a λ-value from that vector, the factorization will be used (otherwise it will be computed).

Further recycling is possible. If the variable `max_factorizations` is set
to a positive value, the object will store that many factorizations
for possible reuse. Every `lin_solve`-call then computes a factorization,
unless a `lin_solve`-call for that `λ` has been computed earlier.
This procedure can at most store `max_factorization` (which can be set `Inf`).

See also: [`FactorizeLinSolver`](@ref), [`create_linsolver`](@ref)

"""
struct FactorizeLinSolverCreator{T_values,T_factor} <: LinSolverCreator
    umfpack_refinements::Int;
    recycled_factorizations::Dict{T_values,T_factor};
    max_factorizations::Int
    function FactorizeLinSolverCreator(;umfpack_refinements::Int=2,
                                       max_factorizations=0,
                                       nep=nothing,
                                       precomp_values=[]
                                       )

        if (precomp_values isa Number) # Helper: Allow values to be a scalar. Put in a vector
            precomp_values=[precomp_values];
        end

        if (size(precomp_values,1)>0 && nep==nothing)
            error("When you want to precompute factorizations you need to supply the keyword argument `nep`");
        end

        # Compute all the factorizations
        precomp_factorizations=map(s-> factorize(compute_Mder(nep,s)), precomp_values);


        # Put them in a dict
        T_from=eltype(precomp_values);
        local T_to,T_from
        if (size(precomp_values,1)>0)
            T_to=eltype(precomp_factorizations);
        else
            T_to=Any;
            T_from=Number;
        end

        dict=Dict{T_from,T_to}();
        for i=1:size(precomp_values,1)
            dict[precomp_values[i]]=precomp_factorizations[i];
        end

        return new{T_from,T_to}(umfpack_refinements,dict,max_factorizations)

    end
end


function create_linsolver(creator::FactorizeLinSolverCreator,nep,λ)
    # Let's see if we find it in recycled_factorizations
    if (λ in keys(creator.recycled_factorizations))
        Afact=creator.recycled_factorizations[λ];
        return FactorizeLinSolver(Afact,creator.umfpack_refinements);
    else
        solver=FactorizeLinSolver(nep,λ,creator.umfpack_refinements);
        if  (length(keys(creator.recycled_factorizations))
             < creator.max_factorizations )
            # Let's save the factorization
            creator.recycled_factorizations[λ]=solver.Afact;
        end
        return solver;
    end

end

struct GMRESLinSolverCreator{T} <: LinSolverCreator where {T}
    kwargs::T
end

"""
    GMRESLinSolverCreator(;kwargs...)

This is the creator for the GMRES-method. Instantiate this
object if you want to use GMRES as your linear system solver.
The `kwargs` are stored and used as keyword arguments in the
call to gmres. See list
of keyword in the [IterativeSolvers.jl manual](https://juliamath.github.io/IterativeSolvers.jl/dev/linear_systems/gmres/).

"""
function GMRESLinSolverCreator(;kwargs...)
    return GMRESLinSolverCreator{typeof(kwargs)}(kwargs)
end


function create_linsolver(creator::GMRESLinSolverCreator,nep, λ)
    return GMRESLinSolver{typeof(λ)}(nep, λ, creator.kwargs)
end


"""
    DeflatedNEPLinSolverCreator(orglinsolvercreator)

This is the creator for case of a deflated NEP.
The argument `orglinsolvercreator` is the `LinSolverCreator` for the original NEP.
The extended linear system

```math
[M, U; X^T, 0][v1; v2] = [b1; b2]
```

is solved with a Schur complement strategy, recycling
the linear solver of the original NEP. Hence, pre-computed entities such as, e.g.,
factorizations and preconditioners can be reused.


NB1: The implementation assumes minimality index = 1.
NB2: The Schur complement is explicitly formed. Hence it is only efficient for a
few deflated eigenvalues.

See also: [`DeflatedNEPLinSolver`](@ref), [`create_linsolver`](@ref),
[`deflate_eigpair`](@ref)

# References
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
"""
struct DeflatedNEPLinSolverCreator{T_origcreator} <: LinSolverCreator where {T_origcreator <: LinSolverCreator}
    orglinsolvercreator::T_origcreator
end

function create_linsolver(creator::DeflatedNEPLinSolverCreator, nep::DeflatedNEP, λ)
    orglinsolver = create_linsolver(creator.orglinsolvercreator, nep.orgnep, λ)
    return DeflatedNEPLinSolver(nep, λ, orglinsolver)
end


# For the moment, Factorize is the default behaviour
"""
    DefaultLinSolverCreator

This is the default linear solver if no other is specified (for most methods).
It is a `FactorizeLinSolverCreator`.


See also: [`LinSolver`](@ref), [`create_linsolver`](@ref),
[`lin_solve`](@ref), [`FactorizeLinSolverCreator`](@ref), [`FactorizeLinSolver`](@ref)

"""
DefaultLinSolverCreator = FactorizeLinSolverCreator
