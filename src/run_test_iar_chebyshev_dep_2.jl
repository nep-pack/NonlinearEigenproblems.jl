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

nep=SPMF_NEP([eye(n), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])

# Then it is needed to create a type to access to this function
import NEPSolver.ComputeY0Cheb
import NEPSolver.ComputeY0ChebPEP
import NEPSolver.ComputeY0ChebDEP
import NEPSolver.AbstractPrecomputeData
import NEPSolver.precompute_data
import NEPSolver.compute_y0_cheb


abstract type ComputeY0Cheb_QDEP <: NEPSolver.ComputeY0Cheb end
type PrecomputeData_QDEP <: AbstractPrecomputeData
    precomp_PEP; precomp_DEP; nep_pep; nep_dep
end

function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    # split the problem as PEP+DEP
    # M(λ)=-λ^2 I + A0 + A1 exp(-τ λ) =
    #     = (A0 + λ I - λ^2 I) + (-λ I + A1 exp(-τ λ))

    nep_pep=PEP([eye(T,n,n), eye(T,n,n), -eye(T,n,n)]); # the PEP part is defined as
    nep_dep=DEP([A0,A1],[0.0,1.0]);     # the DEP part is defined as
    precomp_PEP=precompute_data(T,nep_pep,NEPSolver.ComputeY0ChebPEP,a,b,m,γ,σ);
    precomp_DEP=precompute_data(T,nep_dep,NEPSolver.ComputeY0ChebDEP,a,b,m,γ,σ);

    # combine the precomputations
    return PrecomputeData_QDEP(precomp_PEP,precomp_DEP,nep_pep,nep_dep)

end

function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
    return compute_y0_cheb(T,precomp.nep_pep,ComputeY0ChebPEP,x,y,M0inv,precomp.precomp_PEP)+ compute_y0_cheb(T,precomp.nep_dep,ComputeY0ChebDEP,x,y,M0inv,precomp.precomp_DEP)+y*(view(precomp.precomp_DEP.Tc,1:size(x,2)+1));
    # y*Tc subtracted at each call of compute_y0iar_cheb. Therefore since we do two calls, we need to add it back once.
end

#TODO: fix this function
v0=ones(n);

errormeasure=default_errmeasure(nep);   # TODO: understand why one needs to do this

λ2,Q2,err2,V2, H2 = iar_chebyshev(nep,maxit=mm,Neig=10,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0,errmeasure=errormeasure);

λ,Q,err,V,H = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=mm,Neig=10,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0,errmeasure=errormeasure);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end


println(norm(V[:,1:10]-V2[:,1:10]))
println(norm(H[1:10,1:10]-H2[1:10,1:10]))
