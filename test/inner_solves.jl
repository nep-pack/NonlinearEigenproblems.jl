# Run tests for the inner solves

workspace()
using NonlinearEigenproblems: NEPCore, NEPTypes, NEPSolver, Gallery
using Base.Test

dep=nep_gallery("dep0",200);
n=size(dep,1);
@testset "Inner Solves" begin
    nn=norm(compute_Mder(dep,0));
    errmeasure= (λ,v) -> norm(compute_Mlincomb(dep,λ,v))/nn;

    pnep=create_proj_NEP(dep);
    Q,R=qr(randn(n,5));
    set_projectmatrices!(pnep,Q,Q)

    λv,V=inner_solve(NEPSolver.DefaultInnerSolver,pnep,λv=[0 1]);
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    λv,V=inner_solve(NEPSolver.NewtonInnerSolver,pnep,λv=[0 1],V=eye(5,2));
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

end
