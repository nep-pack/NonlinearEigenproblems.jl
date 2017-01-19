#  Preliminary code to link with NLEVP toolbox
# 

workspace()
using MATLAB  # Dependence required for successive linear problems
using NEPCore
using NEPSolver

type NLEVP_NEP <: NEP
    n::Integer
    name::String
    Ai::Array
    function NLEVP_NEP(name)
        @mput name
        @matlab begin
            addpath("/home/jarl/jobb/src/nlevp3")
            Ai,funs=nlevp(name)
        @matlab end
        @mget Ai # fetch the matrices
        this=new(size(Ai[1],1),name,Ai);
    end
end

function compute_Mder(nep::NLEVP_NEP,λ::Number,i::Integer=0)
    if (i==0 || i==1)
        M=zeros(nep.Ai[1]);
        (fv,fpv)=call_current_fun(λ)
        local f
        if (i==0)
            f=fv;
        else
            f=fpv
        end
        for i=1:length(nep.Ai)
            println("value:",i," is ",abs(f[i]))
            M=M+nep.Ai[i]*f[i]
        end
        return M
    else
        error("Higher order derivatives not implemented")
    end
end

function call_current_fun(lambda)
    lambda=Complex64(lambda)  # avoid type problems
    @mput lambda
    @matlab begin
        (f,fp)=funs(lambda)
    @matlab end
    @mget f fp
    return f,fp
end

function compute_Mlincomb(nep::NLEVP_NEP,λ::Number,V;a=ones(size(V,2)))
    return compute_Mlincomb_from_Mder(nep,λ,V,a)
end

nep=NLEVP_NEP("gun")
#println("n=",nep.n)
#println("name=",nep.name)
#println("first:")
#M=compute_Mder(nep,10)
#println(norm(M,1))
#println("second:")
#M=compute_Mder(nep,12,1)
#println(norm(M,1))

λ,v=newton(nep)
#    
#function load_it(lambda)
#    @mput lambda
#    @matlab begin
#        addpath("/home/jarl/jobb/src/nlevp3")
#        nn,funs=nlevp("gun")
#    @matlab end
#    @matlab begin
#        f,fp=funs(10)
#    @matlab end
#    @mget f fp
#
#    println("f=",f)
#    println("fp=",fp)
#    #println("fp=",fp)
#    #println("fpp=",fpp)    
#            
#    return f
#end
#
#mynn=load_it(200)
#1
