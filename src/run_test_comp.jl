workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

pep = nep_gallery("pep0");

E,A = companion(pep);

Dc,Vc = eig(A,E);

d = size(pep.A,1)-1;
n = size(pep,1);

x = Vc[(d-1)*n+1:d*n,1];
V = Vc[:,1]; 
λ = Dc[1,1]

print("Computed eigenvalue is: ")
print(λ)
print("\n")

lin_res = norm((λ*E-A)*V);#Norm of the residual of the linearized eigenvalue
print("Residual norm of the linearized eigenvalue problem: ")
print(lin_res)
print("\n")

print("Residual norm of the original PEP: ")
print(norm(compute_Mlincomb(pep,λ,x)))
print("\n")



