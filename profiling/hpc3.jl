# Script to be run as profiling on HPC setting

using LinearAlgebra
using NonlinearEigenproblems

using GalleryWaveguide

function run_hpc(;preprun=false)
    k=5;  # determines problem size
    nz_nz0_factor=0.1;  # determines preconditioner fineness (as a factor of nz). 0.3 or larger leads to the finest preconditioner.


    primes=[3,5,7,11,13,17, 19]
    #        23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,]
    # Allowed sizes
    szvec=unique(sort(vec(primes*primes')));

    o=ones(Int,size(primes))
    nzvals=kron(kron(kron(primes,primes),primes),primes)

    # Ugly way to determine possible choices for nz0
    nz0vals1=[kron(kron(kron(primes,primes),primes),o)   kron(kron(kron(primes,primes),o),primes) kron(kron(kron(primes,o),primes),primes)  kron(kron(kron(o,primes),primes),primes)]
    #
    nz0vals2=[kron(kron(kron(o,o),primes),primes)  kron(kron(kron(primes,o),o),primes)  kron(kron(kron(primes,primes),o),o)  kron(kron(kron(primes,o),primes),o)  kron(kron(kron(o,primes),primes),o) kron(kron(kron(o,primes),o),primes)]
    #
    nz0vals3=[kron(kron(kron(o,o),o),primes) kron(kron(kron(o,o),primes),o)   kron(kron(kron(o,primes),o),o)   kron(kron(kron(primes,o),o),o)   ]

    nz0vals=[nz0vals1 nz0vals2 nz0vals3]

    II=sortperm(nzvals)
    nzvals=nzvals[II];
    nz0vals=nz0vals[II,:];
    #k1=70; k2=70; # Determines the size


    nz=nzvals[k];
    # Find "closest" nz0
    nz0=Int(nz0vals[k,argmin(abs.(nz0vals[k,:] .- nz*nz_nz0_factor))])


    nx=nz+4
    println("nx=",nx)
    nep=nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = "JARLEBRING", discretization = "FD", neptype = "WEP")

    γ=-3-3.5im
    #γ = -1.3-0.31im
    println("Generating preconditioner nz0=",nz0);
    precond = wep_generate_preconditioner(nep, nz0, γ)

    gmres_kwargs = ((:maxiter,100), (:restart,100), (:log,true), (:Pl,precond), (:tol, 1e-13), (:verbose,false))

    n=size(nep,1);
    println("nep size n=",n);
    v0=ones(n); v0=v0/norm(v0);

    tol=preprun ? 1e-1 : 1e-10;

    λ,v = resinv(ComplexF64,nep,displaylevel=1,λ=γ,v=v0,
                 tol=tol,
                 linsolvercreator=(nep,γ) -> wep_gmres_linsolvercreator(nep, γ, gmres_kwargs))

end
