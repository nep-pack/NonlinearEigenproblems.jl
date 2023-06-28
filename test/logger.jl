# Tests for NEP types functionality

using NonlinearEigenproblems
using Test

@bench @testset "logger (ErrorLogger)" begin
    nep=nep_gallery("dep0")
    logger=ErrorLogger(1,100,0);
    λ1,x1 = newton(nep,logger=logger,v=ones(size(nep,1)),λ=0,tol=eps()*10)
    E=logger.errs[:,1];
    EE=E[(!isnan).(E)];
    # Check "essentially" quad convergence
    @test abs(log10(EE[end-1]^2)-log10(EE[end]))<2
    @test length(EE) > 5 #Has done more than 5 iterations
    @test EE[end] < eps()*10 # Error measure fulfills stopping criteria


end
