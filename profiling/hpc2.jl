# Script to be run as profiling on HPC setting

using LinearAlgebra
using NonlinearEigenproblems

using GalleryWaveguide

import GalleryWaveguide.SchurMatVec

function run_hpc(;preprun=false)
    primes=[3,5,7,11,13,17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,]
    # Allowed sizes
    szvec=unique(sort(vec(primes*primes')));

    k1=70; k2=70; # Determines the size

    nx=szvec[k1];
    nz=szvec[k2];
    nep=nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = "JARLEBRING", discretization = "FD", neptype = "WEP")




    n=size(nep,1);
    println("nep size n=",n);
    λ0=-3-3.5im
    v0=ones(n); v0=v0/norm(v0);
    tol=preprun ? 1e-1 : 1e-10;
    λ,v = resinv(ComplexF64,nep,displaylevel=1,λ=λ0,v=v0,
                 tol=tol)

end
