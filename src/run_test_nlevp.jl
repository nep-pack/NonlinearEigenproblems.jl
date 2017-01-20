#  Example illustrating loading a NEP from the NLEVP toolbox
# 

workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using MATLAB  
using NEPCore
using NEPSolver
using NEPSolver_MSLP
# We have to explicitly specify functions that we want "overload"
import NEPCore.compute_Mder
import NEPCore.size



"""
     NLEVP_NEP represents a NEP in the NLEVP-toolbox
Example usage: nep=NLEVP_NEP("gun")
"""
type NLEVP_NEP <: NEP
    n::Integer
    name::String
    Ai::Array
    function NLEVP_NEP(name;nlevp_path::String="../../nlevp3")
        if (~isfile(joinpath(nlevp_path,"nlevp.m")))
            error("nlevp.m not found. You need to install the Berlin-Manchester collection (http://www.maths.manchester.ac.uk/our-research/research-groups/numerical-analysis-and-scientific-computing/numerical-analysis/software/nlevp/) and specify a nlevp_path.")
        end
        @mput name nlevp_path
        @matlab begin
            addpath(nlevp_path)
            Ai,funs=nlevp(name)
        @matlab end
        @mget Ai # fetch and store the matrices
        this=new(size(Ai[1],1),name,Ai);
    end
end

function compute_Mder(nep::NLEVP_NEP,λ::Number,i::Integer=0)
    if (i==0 || i==1)        
        lambda=Complex{Float64}(λ)  # avoid type conversion problems
        #println("type",typeof(lambda))
#        nep_name::String=nep.name
#        @mput lambda nep_name
#        if (i==0)
#            println(λ)
#            @matlab begin
#                ll=1+0.1i
#                M=nlevp("eval",nep_name,lambda
#            @matlab end
#            @mget M
#            return M
#        else
#            @matlab begin
#                (M,Mp)=nlevp("eval",nep_name,lambda)
#            @matlab end
#            @mget Mp
#            return Mp
#        end
#    return f,fp
        (fv,fpv)=call_current_fun(lambda)
        M=zeros(nep.Ai[1]);
        if (i==0)
            f=fv;
        else
            f=fpv
        end
        for i=1:length(nep.Ai)
            #println("value:",i," is ",abs(f[i]))
            M=M+nep.Ai[i]*f[i]
        end
        return M
    else
        error("Higher order derivatives not implemented")
    end
end

# Return function and derivative if the current matlab function funs
function call_current_fun(lambda)
    lambda=Complex64(lambda)  # avoid type problems
    @mput lambda
    @matlab begin
        (f,fp)=funs(lambda)
    @matlab end
    @mget f fp
    return f,fp
end

function size(nep::NLEVP_NEP,dim=-1)
    if (dim==-1)
        return (nep.n,nep.n)
    else
        return nep.n
    end
end


nep_name="gun"
println("Loading \"",nep_name,"\" from Berlin-Manchester collection")
nep=NLEVP_NEP(nep_name)
println("Loading completed. Size of problem n=",nep.n)
println("Testing some basic operations on the nep")
M=compute_Mder(nep,160^2+3im)
Mp=compute_Mder(nep,160^2+3im,1)
println("It worked.")


# Approximate eigenvalue of gun problem from Liao, Bai, Lee, Ko
λ0=((1.4948 + 0.000021)*1e2)^2

v0=ones(size(nep,1))

println("Running aug newton")
λ,v=aug_newton(nep,λ=λ0,v=v0,
               displaylevel=2,maxit=30,tolerance=1e-6)
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))
println("Running newton")
λ,v=newton(nep,λ=λ0,v=v0,
               displaylevel=2,maxit=30,tolerance=1e-6)
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))

println("Running MSLP")
λ,v=mslp(nep,λ=λ0,
         displaylevel=2,maxit=30,tolerance=1e-4,eigsolver="inv_it")
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))



println("Running one step of MSLP and then aug_newton")
λ1=λ0
v1=copy(v0)
try 
    λ1,v1=mslp(nep,λ=λ0,
               displaylevel=2,maxit=1,tolerance=-1,
               eigsolver="inv_it")
catch e
    λ1=e.λ
    v1=copy(e.v)
end
v1=v1/norm(v1)
#println("sqrt(λ1)=",sqrt(λ1))
λ,v=aug_newton(nep,λ=λ1,v=v1,c=v0,
               displaylevel=2,maxit=30,tolerance=1e-6)

println("Found eigenvalue \sqrt{λ}=",sqrt(λ))
