using NonlinearEigenproblems, LinearAlgebra, Revise

dep=nep_gallery("dep0",200);
n=size(dep,1);

nn=opnorm(compute_Mder(dep,0));
errmeasure= (λ,v) -> norm(compute_Mlincomb(dep,λ,v))/nn;

pnep=create_proj_NEP(dep);

k = 5
Q,R=qr(randn(n,k));
Q = Matrix(Q)
set_projectmatrices!(pnep,Q,Q)

logger=ErrorLogger(1,100,0);
λv,V = inner_solve(DefaultInnerSolver(), ComplexF64, pnep; λv=[0.0,1.0] .+ 0im, neigs=3, V=ones(k, 2), tol=eps()*100, inner_logger=logger)


logger=ErrorLogger(1,100,0);
λv,V = inner_solve(NewtonInnerSolver(), ComplexF64, pnep; λv=[0.0,1.0] .+ 0im, V=ones(k, 2), tol=eps()*100, inner_logger=logger)

λv,V = inner_solve(IARChebInnerSolver(), ComplexF64, pnep; λv=[0,1,2,3] .+ 0.0im)

λv,V = inner_solve(ContourBeynInnerSolver(), ComplexF64, pnep; λv=[0,1] .+ 0.0im, neigs=3)
