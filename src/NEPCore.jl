module NEPCore

export NEP

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


#NEP() = NEP(4,NaN)
#NEP(n,Md) = NEP(n,Md,default_error_measure)


end
