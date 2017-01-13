module NEPCore

export NEP
export NEP_new
export DEP
export size
export NoConvergenceException
export LinSolver
export compute_Mder
export compute_Mlincomb
export compute_MM

import Base.size  # Overload for nonlinear eigenvalue problems
#using Combinatorics; 

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
    error("No Mder implemented")
    return 0;
end

function compute_Mlincomb(nep::NEP_new,λ::Number,V;a=ones(size(V,2)))
    if (@method_concretely_defined(compute_MM,nep))
        println("Using default: compute_MM -> compute_Mlincomb")
        k=size(V,2)
        S=jordan_matrix(k,λ).'
        b=zeros(size(a));
        for i=1:k
            b[i]=a[i]*factorial(i-1)
        end
        V=V*diagm(b)
        z=compute_MM(nep,S,V)*eye(k,1)
        if (false) # activate for debugging (verify that Mder is equivalent)
           z2=zeros(size(nep,1))
           for i=1:size(a,1)
               z2+=compute_Mder(nep,λ,i-1)*(V[:,i]*a[i])
           end
           println("should be zero:",norm(z2-z))
        end
        return z
    elseif (@method_concretely_defined(compute_Mder,nep))
        # Naive Mlincomb
        println("Using default: compute_MM -> compute_Mlincomb")
        z=zeros(size(nep,1))
        for i=1:size(a,1)
            z+=compute_Mder(nep,λ,i-1)*(V[:,i]*a[i])
        end
    else
        error("No procedure to compute Mlincomb")
    end
end

function compute_MM(nep::NEP_new,S,V)
    error("No procedure to compute MM")
end


# Overload size function. Note: All NEPs must have a field: n.
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



"""
    compute_Mder(nep::DEP,λ::Number,i::Integer=0)
 Compute the ith derivative of a DEP
"""
function compute_MM(nep::DEP,S,V)
    Z=-V*S;
    for j=1:size(nep.A,1)
        Z+=nep.A[j]*V*expm(-nep.tauv[j]*S)
    end
    return Z
end
#



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


    function jordan_matrix(n::Integer,λ::Number)
        Z=λ*eye(n)+diagm(ones(n-1),1);
    end
    
end  # End Module
