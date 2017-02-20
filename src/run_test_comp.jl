workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPCore
using LinSolvers
using NEPTypes
using Gallery

pep = nep_gallery("pep0");
d = size(pep.A,1)-1;
n = size(pep,1);

λa=NaN;
xa=NaN;
try
    λa,xa =newton(pep,maxit=30,displaylevel=1);
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λa=e.λ
    xa=e.v
end

print("##### Results for the problem pep0 #####\n")
print("Eigenvalue computed by newton: ",λa,"\n")

Dc,Vc = polyeig(pep,DefaultEigSolver);

ind = 1;
for i=1:d*n
    if(norm(λa-Dc[i]) < 1e-12)
        ind = i
        break
    end
end
print("Closest eigenvalue computed from companion linearization: ",Dc[ind],"\n\n")



#############################################################################################
#############################################################################################

pep = nep_gallery("pep0_sparse_003");
d = size(pep.A,1)-1;
n = size(pep,1);

λa=NaN;
xa=NaN;
try
    λa,xa =newton(pep,maxit=30,λ=-0.75,v=ones(size(pep,1),1),displaylevel=1);
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λa=e.λ
    xa=e.v
end

print("##### Results for the problem pep0_sparse_003 #####\n")
print("Eigenvalue computed by newton: ",λa,"\n")

Dc,Vc = polyeig(pep,SpEigSolver);

ind = 1;
for i=1:d*n
    if(norm(λa-Dc[i]) < 1e-12)
        ind = i
        break
    end
end
print("Closest eigenvalue computed from companion linearization: ",Dc[ind],"\n\n")

#################################################################################
#################################################################################

println("Testing companion linearization with BigFloat");
A0=(Array{BigFloat,2}([1 2; 3 4]));
A1=(Array{BigFloat,2}([1 44; 3 4]));
A2=(Array{BigFloat,2}([1 44; -3 100]));
    
pep3=PEP([A0,A1,A2])
E,A = companion(pep3);
# Power method for testing (since eig does not work for BigFloats)

z=ones(BigFloat,size(A,1));
local evp::BigFloat
local evp_old::BigFloat
local d::BigFloat
TOL=big"1e-50"
d=Inf; evp=Inf
k=0;
while abs(d)>TOL
    k=k+1
    z=z/norm(z);
    z2=E\A*z
    evp_old=evp
    evp=dot(z,z2)
    d=evp-evp_old
    println("Iteration ",k,", λ=",Float64(evp), " Δ = ",Float64(d))
    z=z2
end

println("Solving same problem with res_inv")
λ,v=res_inv(pep3,λ=(BigFloat(evp)+0.1),v=z[1:size(pep3,1)],tolerance=TOL)
#Dc,Vc = compsolve(pep3);


println("Difference btw companion solver and res_inv for bigfloats:",abs(λ-evp))
