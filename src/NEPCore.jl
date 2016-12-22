module NEPCore

export NEP
export NoConvergenceException
export LinSolver

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
            return this.resnorm(λ,v)/norm(Md(λ,0));
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
