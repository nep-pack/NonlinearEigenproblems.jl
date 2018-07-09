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

mm=80;  # number of iterations

nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])

# Then it is needed to create a type to access to this function
import NEPSolver.ComputeY0Cheb
import NEPSolver.AbstractPrecomputeData
abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
type PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_PEP; precomp_DEP
end

# And then introduce a function dispatch for this new type in order to use
# the defined orthogonalization function
import NEPSolver.precompute_data
# function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
#     return PrecomputeData_QDEP()
# end

function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (O I + λ I - λ^2 I) + (-λ I + A0 + A1 exp(-τ λ))

    # the PEP part is defined as
    nep_pep=PEP([zeros(T,n,n), eye(T,n,n), -eye(T,n,n)]);
    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);

    # the DEP part is defined as
    A0=nep.A[2]; A1=nep.A[3];
    nep_dep=DEP([A0,A1],[0,1]);
    precomp_DEP=precompute_data(T,nep_dep,NEPSolver.ComputeY0ChebDEP,a,b,m,γ,σ);

    # combine the precomputations
    return PrecomputeData_QDEP(precomp_PEP,precomp_DEP)

end


import NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (O I + λ I - λ^2 I) + (-λ I + A0 + A1 exp(-τ λ))
    nep_pep=PEP([zeros(T,n,n), eye(T,n,n), -eye(T,n,n)]);
    A0=nep.A[2]; A1=nep.A[3]; τ=1;
    nep_dep=DEP([A0,A1],[0,τ]);

    return compute_y0_cheb(T,nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)+ compute_y0_cheb(T,nep_dep,NEPSolver.ComputeY0ChebDEP,x,y,M0inv,precomp.precomp_DEP)
end

#TODO: fix this function
v0=randn(n);
λ,Q,err = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=mm,Neig=10,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
