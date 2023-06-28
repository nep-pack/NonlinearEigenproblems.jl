# Test of companion linearization

using NonlinearEigenproblems
using Test,Random
using LinearAlgebra

@bench @testset "Companion Linearization" begin

    #####################
    # Dense matrix test #
    #####################
    pep = nep_gallery("pep0");
    deg = size(get_fv(pep),1) -1;
    n = size(pep,1);
    tolerance = 1e-14;

    λa,xa = newton(pep, maxit=30, logger=0, tol = tolerance,v=ones(size(pep,1)));
    Dc,Vc = polyeig(pep, DefaultEigSolver);

    ind = argmin(abs.(λa .- Dc)./abs(λa));
    @test (abs(λa-Dc[ind])/abs(λa) < tolerance*10)


    ######################
    # Sparse matrix test #
    ######################
    pep = nep_gallery("pep0_sparse", 200, 0.03);
    deg = size(get_fv(pep),1) -1;
    n = size(pep,1);
    tolerance = 1e-14;

    λa,xa =newton(pep, maxit=30, λ=-0.75, v=ones(n), logger=0, tol = tolerance*1e-3);
    Dc,Vc = polyeig(pep,DefaultEigSolver);

    ind = argmin(abs.(λa .- Dc)./abs(λa));
    @test (abs(λa-Dc[ind])/abs(λa) < tolerance*4000)


    ########################
    # BigFloat matrix test #
    ########################
    A0 = (Array{BigFloat,2}([1 2; 3 4]));
    A1 = (Array{BigFloat,2}([1 44; 3 4]));
    A2 = (Array{BigFloat,2}([1 44; -3 100]));
    pep = PEP([A0,A1,A2])

    # Power method for testing (since eig does not work for BigFloats)
    E,A = companion(pep);
    z=ones(BigFloat,size(A,1));
    local evp::BigFloat
    local evp_old::BigFloat
    local d::BigFloat
    tolerance = big"1e-50"
    d=Inf; evp=Inf
    k=0;
    while abs(d) > tolerance
        k=k+1
        z=z/norm(z);
        z2=E\A*z
        evp_old=evp
        evp=dot(z,z2)
        d=evp-evp_old
        @debug "Iteration $k, λ=$(Float64(evp)), Δ = $(Float64(d))"
        z=z2
    end
    @info "Solving same problem with resinv"
    λ,v = resinv(BigFloat, pep, λ = (BigFloat(evp)+0.1), v = z[1:size(pep,1)], logger=0, tol = tolerance, errmeasure=ResidualErrmeasure(pep))

    @test (abs(λ-evp)/abs(λ) < tolerance*10)

end
