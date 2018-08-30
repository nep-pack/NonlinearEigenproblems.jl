module NEPCore_old

using LinearAlgebra

export NEP
export size
export NoConvergenceException
export LinSolver

# Core interfaces
export compute_Mder
export compute_Mlincomb
export compute_MM

# NEP-functions
export compute_resnorm

# Helper functions
export compute_Mlincomb_from_MM

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
        Z = λ*eye(n) + diagm(1 => ones(n-1))
    end

end  # End Module
