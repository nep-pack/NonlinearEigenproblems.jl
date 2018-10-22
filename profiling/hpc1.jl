# Script to be run as profiling on HPC setting

using LinearAlgebra
using NonlinearEigenproblems


function run_hpc(;preprun=false)
    nep=nep_gallery("nlevp_native_gun");
    n=size(nep,1)
    S=150^2*[1.0 0; 0 1]; V=[[1 0; 0 1]; zeros(n-2,2)];
    println("blocknewton");
    tol=preprun ? 1e-1 : 1e-10;
    (Z,X)=blocknewton(nep,S=S,X=V,displaylevel=1,armijo_factor=0.5,maxit=20,tol=tol)
end
