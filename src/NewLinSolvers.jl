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


struct FactorizeLinSolverCreator{T_values,T_factor} <: LinSolverCreator
    umfpack_refinements::Int;
    recycled_factorizations::Dict{T_values,T_factor};
    max_factorizations::Int
    function FactorizeLinSolverCreator(;umfpack_refinements::Int=1,
                                       max_factorizations=0,
                                       nep=nothing,
                                       precomp_values=[]
                                       )

        if (size(precomp_values,1)>0 && nep==nothing)
            error("When you want to precompute factorizations you need to supply the keyword argument `nep`");
        end

        # Compute all the factorizations
        local precomp_factorizations

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
# For the moment, Factorize is the default behaviour
DefaultLinSolverCreator = FactorizeLinSolverCreator


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
function GMRESLinSolverCreator(;kwargs...)
    return GMRESLinSolverCreator{typeof(kwargs)}(kwargs)
end


function create_linsolver(creator::GMRESLinSolverCreator,nep::NEP, λ)
    return GMRESLinSolver{typeof(λ)}(nep, λ, creator.kwargs)
end
