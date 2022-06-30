using NonlinearEigenproblems
using Test
using LinearAlgebra
using SpecialFunctions
using SparseArrays


@testset "compute_types" begin
    A0=(1:5)*(1:5)'
    A1=(1:5)*(3:7)'+I;
    nep=DEP([A0,A1]);
    x=ones(5);
    s=compute_rf(ComplexF64,nep,x,
                 NonlinearEigenproblems.NEPSolver.ScalarNewtonInnerSolver(),λ=-0.5+3im)[1];

    @test abs(x'*compute_Mlincomb(nep,s,x)) < eps()*100;

    y=zeros(5); y[1]=3;
    s=compute_rf(ComplexF64,nep,x,
                 NonlinearEigenproblems.NEPSolver.ScalarNewtonInnerSolver(),λ=-0.5+3im,
                 y=y)[1];

    @test abs(y'*compute_Mlincomb(nep,s,x)) < eps()*100;
end
