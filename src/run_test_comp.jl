workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPSolver_MSLP 	
using NEPCore
using NEPTypes
using Gallery

pep = nep_gallery("pep0");

E,A = companion(pep);

Dc,Vc = eig(A,E);
Dc

d = size(pep.A,1)-1;
n = size(pep,1);

x = Vc[(d-1)*n+1:d*n,1];
V = Vc[:,1]; 
λ = Dc[1,1]

λa=NaN;
xa=NaN;
try
    λa,xa =newton(pep,maxit=30,displaylevel=0);
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λa=e.λ
    xa=e.v
end

#println(norm(compute_Mlincomb(pep,λa,xa)))
print("Eigenvalue computed by newton: ",λa,"\n\n")

ind = 1;
for i=1:400
	if(norm(λa-Dc[i]) < 1e-12)
		ind = i
	end
end

print("Closest eigenvalue computed from companion linearization: ")
print(Dc[ind])
	



