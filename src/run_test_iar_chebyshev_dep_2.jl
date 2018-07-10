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


#n=100; A0=rand(n,n); A1=rand(n,n);

mm=80;  # number of iterations

A1=zeros(A1);
nep=SPMF_NEP([eye(n), A0],[λ->-λ^2,λ->eye(λ),λ->expm(-λ])

# Then it is needed to create a type to access to this function
import NEPSolver.ComputeY0Cheb
import NEPSolver.AbstractPrecomputeData
import NEPSolver.precompute_data

abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
import NEPSolver.precompute_data
type PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_PEP; precomp_DEP; nep_pep; nep_dep
end

function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (O I + λ I - λ^2 I) + (-λ I + A0 + A1 exp(-τ λ))


    nep_pep=PEP([A0, eye(T,n,n), -eye(T,n,n)]); # the PEP part is defined as
    nep_dep=DEP([A1],[1.0]);     # the DEP part is defined as


    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);

    precomp_DEP=precompute_data(T,nep_dep,NEPSolver.ComputeY0ChebDEP,a,b,m,γ,σ);

    # combine the precomputations
    return PrecomputeData_QDEP(precomp_PEP,precomp_DEP,nep_pep,nep_dep)

end


import NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    return compute_y0_cheb(T,precomp.nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)+ compute_y0_cheb(T,precomp.nep_dep,NEPSolver.ComputeY0ChebDEP,x,y,M0inv,precomp.precomp_DEP)
end

#TODO: fix this function
v0=ones(n);

λ2,Q2,err2,V2 = iar_chebyshev(nep,maxit=mm,Neig=20,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0);

λ,Q,err,V = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=mm,Neig=20,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
