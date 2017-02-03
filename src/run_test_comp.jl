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


#This residual is of the order 10-12
lin_res = norm((λ*E-A)*V)#Norm of the residual of the linearized eigenvalue

#Gives a very large residual

#Mλ = pep.A[1];
#for i=2:d+1
#    Mλ = Mλ + pep.A[i]*λ^(i-1);
#end


#nlin_res = norm(Mλ*x);#Norm of the residual of the original PEP
