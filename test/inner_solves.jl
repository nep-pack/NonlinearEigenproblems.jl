# Run tests for the inner solves

using NonlinearEigenproblems
using Test,Random
using LinearAlgebra

Random.seed!(0);
@bench @testset "Inner Solves" begin
    dep=nep_gallery("dep0",200);
    n=size(dep,1);
    Random.seed!(0);

    nn=opnorm(compute_Mder(dep,0));
    errmeasure= (λ,v) -> norm(compute_Mlincomb(dep,λ,v))/nn;

    pnep=create_proj_NEP(dep);

    k = 5
    temp_pep=nep_gallery("pep0",200); Q,R=qr(temp_pep.A[1][:,8:12]);  Q = Matrix(Q) # Ugly but gives stability over versions. Stable "random" projection space Q.
    set_projectmatrices!(pnep,Q,Q)

    logger=ErrorLogger(50,100,0);
    λv,V = inner_solve(DefaultInnerSolver(), ComplexF64, pnep; σ=0.0+0im, neigs=3, tol=eps()*100, inner_logger=logger)
    verify_lambdas(3, pnep, λv, V,eps()*10000)
    E=logger.errs[:,1]; E=E[(!isnan).(E)];

    @test length(E) > 5 #Has done more than 5 iterations
    @test E[end] <eps()*10000 # Error measure fulfills stopping criteria

    logger=ErrorLogger(1,100,0);
    λv,V = inner_solve(NewtonInnerSolver(), ComplexF64, pnep; λv=[0.0,1.0] .+ 0im, V=ones(k, 2), tol=eps()*100, inner_logger=logger)
    verify_lambdas(2, pnep, λv, V,eps()*10000)
    E=logger.errs[:,1]; E=E[(!isnan).(E)];
    @test length(E) > 5 #Has done more than 5 iterations
    @test E[end] <eps()*10000 # Error measure fulfills stopping criteria

    #λv,V=inner_solve(SGIterInnerSolver(),pnep,λv=[0.0],j=1);
    #verify_lambdas(2, pnep, λv, V,eps()*10000)

    λv,V = inner_solve(IARChebInnerSolver(), ComplexF64, pnep; λv=[0,1,2,3] .+ 0.0im)
    verify_lambdas(10, pnep, λv, V, sqrt(eps()))

    # TODO: this test results in a "Rank drop" warning, and a third eigenvalue that's not converged
    λv,V = inner_solve(ContourBeynInnerSolver(), ComplexF64, pnep; λv=[0,1] .+ 0.0im, neigs=3)
    # @test minimum(svdvals(compute_Mder(pnep,λv[1]))) < eps()*100
    verify_lambdas(2, pnep, λv[1:2], V[:,1:2], sqrt(eps()))

    λv,V = inner_solve(NleigsInnerSolver(), ComplexF64, pnep; λv=[0,1,2,3] .+ 0.0im)
    @test size(λv,1) >= 5;
end
