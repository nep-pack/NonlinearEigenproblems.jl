module NEPCore

export NEP
export NoConvergenceException

type NEP
    n::Int
    Md
    resnorm
    relresnorm
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
#
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



end
