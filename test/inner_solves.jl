# Run tests for the inner solves

workspace()
using NonlinearEigenproblems: NEPCore, NEPTypes, NEPSolver, Gallery
using Base.Test

#import NEPSolver.inner_solve;
#include("../src/inner_solver.jl");

dep=nep_gallery("dep0",200);
n=size(dep,1);
@testset "Inner Solves" begin
    nn=norm(compute_Mder(dep,0));
    errmeasure= (λ,v) -> norm(compute_Mlincomb(dep,λ,v))/nn;

    pnep=create_proj_NEP(dep);
    Q,R=qr(randn(n,5));
    set_projectmatrices!(pnep,Q,Q)

    λv,V=inner_solve(NEPSolver.DefaultInnerSolver,Complex128,pnep,λv=[0.0,1.0]+0im,Neig=3);
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    λv,V=inner_solve(NEPSolver.NewtonInnerSolver,Complex128,pnep,λv=[0.0,1.0]+0im,V=eye(5,2),tol=eps()*100);
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    #λv,V=inner_solve(NEPSolver.SGIterInnerSolver,pnep,λv=[0.0],j=1);
    #@test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    λv,V=inner_solve(NEPSolver.IARChebInnerSolver,Complex128,pnep,λv=[0,1,2,3]+0.0im);
    nn=norm(compute_Mlincomb(pnep,λv[1],V[:,1]));
    @test nn < eps()*100

    λv,V=inner_solve(NEPSolver.ContourBeynInnerSolver,Complex128,pnep,λv=[0,1]+0.0im,Neig=3);
    nn=minimum(svdvals(compute_Mder(pnep,λv[1])))
    @test nn < eps()*100

end
