# Run tests for the inner solves

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems
using Test
using LinearAlgebra

@bench @testset "Inner Solves" begin
    dep=nep_gallery("dep0",200);
    n=size(dep,1);

    nn=opnorm(compute_Mder(dep,0));
    errmeasure= (λ,v) -> norm(compute_Mlincomb(dep,λ,v))/nn;

    pnep=create_proj_NEP(dep);

    k = 5
    Q,R=qr(randn(n,k));
    Q = Matrix(Q)
    set_projectmatrices!(pnep,Q,Q)

    λv,V = inner_solve(NEPSolver.DefaultInnerSolver, ComplexF64, pnep; λv=[0.0,1.0] .+ 0im, Neig=3, V=ones(k, 2), tol=eps()*100)
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    λv,V = inner_solve(NEPSolver.NewtonInnerSolver, ComplexF64, pnep; λv=[0.0,1.0] .+ 0im, V=ones(k, 2), tol=eps()*100)
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    #λv,V=inner_solve(NEPSolver.SGIterInnerSolver,pnep,λv=[0.0],j=1);
    #@test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    λv,V = inner_solve(NEPSolver.IARChebInnerSolver, ComplexF64, pnep; λv=[0,1,2,3] .+ 0.0im)
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100

    λv,V = inner_solve(NEPSolver.ContourBeynInnerSolver, ComplexF64, pnep; λv=[0,1] .+ 0.0im, Neig=3)
    # @test minimum(svdvals(compute_Mder(pnep,λv[1]))) < eps()*100
    @test norm(compute_Mlincomb(pnep,λv[1],V[:,1])) < eps()*100
end
