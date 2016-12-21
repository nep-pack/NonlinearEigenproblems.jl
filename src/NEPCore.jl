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
    # 
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
            # One step of scalar Newton's method
            λ=λ0;
            for k=1:10
                M=this.Md(λ,0);
                Md=this.Md(λ,1);
                Δλ=-dot(y,M*x)/dot(y,Md*x);
                λ=λ+Δλ
            end
            return λ
        end
        #

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
end
