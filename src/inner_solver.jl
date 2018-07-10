# Helper functions for methods based on inner-outer iterations

abstract type InnerSolver end;

abstract type NewtonInnerSolver end;
abstract type PolyeigInnerSolver end;
abstract type DefaultInnerSolver end;

function inner_solve(TT::Type{T},nep::Proj_NEP;kwargs...) where T<:NewtonInnerSolver

    kvargsdict=Dict(kwargs);
    λv=kv
    argsdict[:λv];
    V =kvargsdict[:V];
    tol =kvargsdict[:tol];
    for k=1:size(λv,1)
        try
            v0=V[:,k]; # Starting vector for projected problem
            projerrmeasure=(λ,v) -> norm(compute_Mlincomb(nep,λ,v))/norm(compute_Mder(nep,λ));
            # Compute a solution to projected problem with Newton's method
            λ1,vproj=newton(nep,displaylevel=0,λ=λv[k],
                            v=v0,maxit=50,tol=tol/10,
                            errmeasure=projerrmeasure);            
            V[:,k]=vproj;
            λv[k]=λ1;
        catch e
            if (isa(e, NoConvergenceException))
                λ1=λv[k];
                vproj=V[:,k];
            else
                rethrow(e)
            end
        end

        
    end
    return λv,V
end

function inner_solve(TT::Type{T},nep::Proj_NEP;kwargs...) where T<:PolyeigInnerSolver
    if (typeof(nep.orgnep)!=PEP)
        error("Wrong type");
    end
    pep=PEP(get_Av(nep.nep_proj))
    return polyeig(pep);
end

    
