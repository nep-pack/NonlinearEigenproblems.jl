workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPSolver_MSLP 	
using NEPCore
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

E,A = companion(pep);

s = LinEigSolver()
Dc,Vc = s.solve(A,B=E,λ_t=0,nev=d*n);

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


E,A = companion(pep);

Dc,Vc = s.solve(A,B=E,λ_t=1,nev=d*n);

ind = 1;
for i=1:d*n
    if(norm(λa-Dc[i]) < 1e-12)
        ind = i
        break
    end
end
print("Closest eigenvalue computed from companion linearization: ",Dc[ind],"\n\n")




