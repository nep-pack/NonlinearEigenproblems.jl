workspace(); push!(LOAD_PATH, pwd()); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver
#import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");

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


n=100; A0=rand(n,n); A1=rand(n,n);
#A0=eye(A0); A1=eye(A1);

mm=80;  # number of iterations

#A1=zeros(A1);
nep=SPMF_NEP([A0, A1],[λ->eye(λ),λ->expm(-λ)])

# Then it is needed to create a type to access to this function
import NEPSolver.ComputeY0Cheb
import NEPSolver.AbstractPrecomputeData
import NEPSolver.precompute_data

abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
import NEPSolver.precompute_data
type PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_DEP; nep_dep
end

function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    nep_dep=DEP([A0,A1],[0.0,1.0]);     # the DEP part is defined as
    precomp_DEP=precompute_data(T,nep_dep,NEPSolver.ComputeY0ChebDEP,a,b,m,γ,σ);
    return PrecomputeData_QDEP(precomp_DEP,nep_dep)
end

T=Complex128; precomp=precompute_data(T,nep,ComputeY0Cheb_QDEP,-1,1,mm,1,0)

import NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    return compute_y0_cheb(T,precomp.nep_dep,NEPSolver.ComputeY0ChebDEP,x,y,M0inv,precomp.precomp_DEP)
end

#TODO: fix this function
v0=ones(n);

errormeasure=default_errmeasure(precomp.nep_dep);
#errormeasure=default_errmeasure(nep);  # TODO: understand why this doesn't work

λ2,Q2,err2,V2,H2 = iar_chebyshev(precomp.nep_dep,maxit=mm,Neig=mm,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0,errmeasure=errormeasure);

λ,Q,err,V,H = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=mm,Neig=mm,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0,errmeasure=errormeasure);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

norm(V-V2)
norm(H-H2)
