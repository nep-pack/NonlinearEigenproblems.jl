workspace(); push!(LOAD_PATH, pwd()); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver
#import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");

#TODO: add this to the gallery
n=100;

mm=80;  # number of iterations

A0=rand(n,n); A1=rand(n,n); A2=rand(n,n);

nep=SPMF_NEP([A0, A1, A2],[λ->eye(λ),λ->λ,λ->λ^2])

# Then it is needed to create a type to access to this function
import NEPSolver.ComputeY0Cheb
import NEPSolver.AbstractPrecomputeData
import NEPSolver.precompute_data

abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
import NEPSolver.precompute_data
type PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_PEP; nep_pep
end


function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (O I + λ I - λ^2 I) + (-λ I + A0 + A1 exp(-τ λ))

    # the PEP part is defined as
    nep_pep=PEP([A0,A1,A2])
    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);

    # combine the precomputations
    return PrecomputeData_QDEP(precomp_PEP,nep_pep)

end


import NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    return compute_y0_cheb(T,precomp.nep_pep,NEPSolver.ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)
end

#TODO: fix this function
v0=randn(n);
λ,Q,err = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=mm,Neig=10,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
