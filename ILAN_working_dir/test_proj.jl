using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise

dep=nep_gallery("dep0",200);
n=size(dep,1);

nn=opnorm(compute_Mder(dep,0));
errmeasure= (λ,v) -> norm(compute_Mlincomb(dep,λ,v))/nn;

pnep=create_proj_NEP(dep);

k = 5
Q,R=qr(randn(n,k));
Q = Matrix(Q)
set_projectmatrices!(pnep,Q,Q)

λv,V = inner_solve(IARInnerSolver(maxit=150), ComplexF64, pnep; neigs=Inf)
λv
