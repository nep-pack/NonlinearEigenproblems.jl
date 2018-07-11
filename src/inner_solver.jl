# Helper functions for methods based on inner-outer iterations

export inner_solve;
export InnerSolver;

abstract type InnerSolver end;
abstract type NewtonInnerSolver <: InnerSolver end;
abstract type PolyeigInnerSolver <: InnerSolver end;
abstract type DefaultInnerSolver <: InnerSolver end;
abstract type IARInnerSolver <: InnerSolver end;
abstract type IARChebInnerSolver <: InnerSolver end;
abstract type SGIterInnerSolver <: InnerSolver end;

"""
  inner_solve(T,nep;kwargs...)

Solves the projected linear problem with solver specied with T. This is to be used
as an inner solver in an inner-outer iteration. T specifies which method
to use. The most common choice is `DefaultInnersolver`. The function returns
`(λv,V)` where `λv` is an array of eigenvalues and `V` a matrix with corresponding
 vectors.
"""
function inner_solve(TT::Type{DefaultInnerSolver},nep::NEPTypes.Proj_NEP;kwargs...)
    if (typeof(nep.orgnep)==NEPTypes.PEP)
        return inner_solve(PolyeigInnerSolver,nep;kwargs...);
    elseif (typeof(nep.orgnep)==NEPTypes.DEP)
        # Should be Cheb IAR
        return inner_solve(IARChebInnerSolver,nep;kwargs...);
    elseif (typeof(nep.orgnep)==NEPTypes.SPMF_NEP) # Default to IAR for SPMF
        return inner_solve(IARInnerSolver,nep;kwargs...);
    else
        return inner_solve(NewtonInnerSolver,nep;kwargs...);        
    end
end


function inner_solve(TT::Type{NewtonInnerSolver},nep::NEPTypes.Proj_NEP;kwargs...) 

    kvargsdict=Dict(kwargs);
    λv=kvargsdict[:λv];
    V =kvargsdict[:V];
    tol =kvargsdict[:tol];
    for k=1:size(λv,1)
        try
            v0=V[:,k]; # Starting vector for projected problem
            projerrmeasure=(λ,v) -> norm(compute_Mlincomb(nep,λ,v))/norm(compute_Mder(nep,λ));
            # Compute a solution to projected problem with Newton's method
            λ1,vproj=augnewton(nep,displaylevel=0,λ=λv[k],
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

function inner_solve(TT::Type{PolyeigInnerSolver},nep::NEPTypes.Proj_NEP;kwargs...) 
    if (typeof(nep.orgnep)!=NEPTypes.PEP)
        error("Wrong type");
    end
    pep=NEPTypes.PEP(NEPTypes.get_Av(nep.nep_proj))
    return polyeig(pep);
end


function inner_solve(TT::Type{IARInnerSolver},nep::NEPTypes.Proj_NEP;λv=[],kwargs...) 
    m=size(λv,1);
    σ=mean(λv);
    try 
        λ,V=iar(nep,σ=σ,Neig=m+3,tol=1e-13,maxit=50);
        return λ,V
    catch e
        if (isa(e, NoConvergenceException))
            # If we fail to find the wanted number of eigvals, we can still
            # return the ones we found
            return e.λ,e.v
        else
            rethrow(e)
        end
    end
end

function inner_solve(TT::Type{IARChebInnerSolver},nep::NEPTypes.Proj_NEP;λv=[],kwargs...) 
    m=size(λv,1);
    σ=mean(λv);
    if (typeof(nep.orgnep)==NEPTypes.DEP)
        AA=get_Av(nep);
        BB=Array{eltype(AA),1}(size(AA,1)-1);
        for i=1:(size(AA,1)-1)
            BB[i]=AA[1]\AA[1+i];
        end
        println("BB=",BB);;
        println("tauv=",nep.orgnep.tauv);
        println("m=",m);;
        println("σ=",σ);;        
        nep=DEP(BB,nep.orgnep.tauv)
    end

    println("typeof(nep)=",typeof(nep));
    try 
        λ,V=iar_chebyshev(nep,σ=σ,Neig=m+3,tol=1e-13,maxit=50);
        return λ,V
    catch e
        if (isa(e, NoConvergenceException))
            # If we fail to find the wanted number of eigvals, we can still
            # return the ones we found
            return e.λ,e.v
        else
            rethrow(e)
        end
    end
end


function inner_solve(TT::Type{SGIterInnerSolver},nep::NEPTypes.Proj_NEP;λv=[0],j=0,kwargs...) 
    m=size(λv,1);
    σ=mean(λv);
    λ,V=sgiter(nep,j)
    return λ,V
#    try 
#        λ,V=sgiter(nep,j)
#        return λ,V
#    catch e
#        if (isa(e, NoConvergenceException))
#            # If we fail to find the wanted number of eigvals, we can still
#            # return the ones we found
#            return e.λ,e.v
#        else
#            rethrow(e)
#        end
#    end
end

