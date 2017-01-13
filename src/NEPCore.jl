module NEPCore

export NEP
export NEP_new
export DEP
export size
export NoConvergenceException
export LinSolver
export compute_Mder

import Base.size  # Overload for nonlinear eigenvalue problems

"""
Determines if a method is defined for the concrete class.
False will be returned if methodname is only defined for the
abstract superclass of s.
"""
macro method_concretely_defined(methodname,s)
    return :(length(methodswith(typeof($s), $methodname, false))>0)
end


abstract NEP_new;


############################################
# Default NEP functions
#
"""
    compute_Mder(nep::NEP_new,λ::Number,i::Integer=0)    
 Computes the ith derivative of NEP evaluated in λ\\
 Usage:\\
   compute_Mder(nep,λ)  # Evaluate NEP in λ
"""
function compute_Mder(nep::NEP_new,λ::Number,i::Integer=0)
    return 3
end

# Overload size function. Note: All NEPs must have a field n.
function size(nep::NEP_new,dim=-1)
    if (dim==-1)
        return (nep.n,nep.n)
    else
        return nep.n
    end
end



"""
    Delay eigenvalue problem
"""
type DEP <: NEP_new
    n::Integer
    A   # A can be a full or sparse matrix
    tauv::Array{Float64,1}
    function DEP(AA,tauv=[0,1.0])
        n=size(AA[1],1)
        this=new(n,AA,tauv);
        return this;
    end
end
"""
    compute_Mder(nep::DEP,λ::Number,i::Integer=0)
 Compute the ith derivative of a DEP
"""
function compute_Mder(nep::DEP,λ::Number,i::Integer=0)
    local M,I;
    if issparse(nep.A[1])
        M=spzeros(nep.n,nep.n)
        I=speye(nep.n,nep.n)
    else
        M=zeros(nep.n,nep.n)
        I=eye(nep.n,nep.n)        
    end
    if i==0; M=-λ*I;  end
    if i==1; M=-I; end
    for j=1:size(nep.A,1)
         M+=nep.A[j]*(exp(-nep.tauv[j]*λ)*(-nep.tauv[j])^i)
    end
    return M
end





### Old stuff

type NEP
    n::Int
    Md
    resnorm
    relresnorm


    # Rayleigh functional:
    # Compute λ s.g. y^TM(λ)x=0
    rf::Function

    # Linear combination of derivatives
    # compute a[1]*M^{(0)}(s)V[:,1]+...+a[p]*M^{(p-1)}(s)V[:,p]
    Mlincomb::Function
    
    function NEP(n,Md)
        this=new()
        this.n=n
        this.Md=Md

        this.resnorm=function (λ,v)
            return norm(this.Md(λ,0)*v)
        end

        this.relresnorm=function (λ,v)
             # Giampaolo: I changed to 1-norm in order to expand 
             # to sparse matrices
            return this.resnorm(λ,v)/norm(Md(λ,0),1);
        end

        this.rf=function(x; y=x, target=0, λ0=target)
            # Ten steps of scalar Newton's method
            λ=λ0;
            for k=1:10
                Δλ=-dot(y,this.Md(λ,0)*x)/dot(y,this.Md(λ,1)*x);
                λ=λ+Δλ
            end
            return λ
        end

        this.Mlincomb=function (λ,V;a=ones(size(V,2)))
            z=zeros(this.n)
            for k=1:size(V,2)
                z+=this.Md(λ,k-1)*V[:,k]*a[k]
            end
            return z
        end
        return this
            
    end
end


# In case an iterative method does not converge
type NoConvergenceException <: Exception
    λ
    v
    errmeasure
    msg
end


# To slove linear systems
type LinSolver
    solve::Function
    Afact
    function LinSolver(A)
        this=new()
        this.Afact=factorize(A)
        this.solve=function foo(x;tol=eps())
            return this.Afact\x
        end
        return this;
    end   
end


end  # End Module
