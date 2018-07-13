#workspace(); push!(LOAD_PATH, pwd()); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver
import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");

n=100; A0=rand(n,n); A1=rand(n,n);

mm=80;  # number of iterations

nep=SPMF_NEP([eye(n), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])
errormeasure=default_errmeasure(nep);   # TODO: understand why one needs to do this

λ,Q,err,V, H = iar_chebyshev(nep,maxit=mm,Neig=20,σ=0.0,γ=1,displaylevel=0,check_error_every=1,errmeasure=errormeasure,v=ones(n));

λ2,Q2,err2,V2, H2 = iar_chebyshev(nep,maxit=mm,Neig=20,σ=0.0,γ=1,displaylevel=0,check_error_every=1,errmeasure=errormeasure,compute_y0_method=ComputeY0Cheb,v=ones(n));


norm(V[:,1:10]-V2[:,1:10])
