#workspace(); push!(LOAD_PATH, pwd()); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver
import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl"); import NEPSolver.iar; include("../src/method_iar.jl"); import NEPCore.compute_Mlincomb_from_MM!; include("../src/NEPCore.jl");

#TODO: add this to the gallery
n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];


#n=100; A0=rand(n,n); A1=rand(n,n);

mm=80;  # number of iterations

A1=zeros(A1);
# something strange happen here
nep=SPMF_NEP([eye(n), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])
#nep=SPMF_NEP([eye(n), A0],[λ->-λ^2,λ->eye(λ)])


λ,Q,err,V = iar_chebyshev(nep,maxit=mm,Neig=8,σ=0.0,γ=1,displaylevel=1,check_error_every=1);

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
