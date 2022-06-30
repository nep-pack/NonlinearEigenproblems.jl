# Tests for core functionality

using NonlinearEigenproblems
using LinearAlgebra
using Test,Random

struct TestNEP <: NEP
   dummy
end

@bench @testset "Core" begin

nep=nep_gallery("dep0");
n=size(nep,1);

## Mlincomb_tests
Random.seed!(0);
V=complex(randn(n,3));
λ=0.3+1im;
T=typeof(λ);
z1=compute_Mlincomb(nep,λ,V)
z2=compute_Mlincomb(nep,λ,V,[T(1),T(1),T(1)])
@test z1==z2
z3=compute_Mder(nep,λ,0)*V[:,1]+ compute_Mder(nep,λ,1)*V[:,2]+ compute_Mder(nep,λ,2)*V[:,3]
@test z1 ≈ z3

## Test Mlincomb starting counting with kth derivative
z2=compute_Mlincomb(nep,λ,V,[T(0),T(0),T(1)])
z3=compute_Mder(nep,λ,2)*V[:,3]
@test z2≈z3
z4=compute_Mlincomb(nep,λ,V[:,3],[T(1)], 2);
@test z3≈z4
#

## Simple coverage. Test that methods throw error if not defined
my_test_NEP = TestNEP([1 2; 3 4])
@test_throws ErrorException compute_Mder(my_test_NEP, 1+1im, 2)
#@test_throws MethodError compute_Mlincomb(my_test_NEP , 1+1im, [1 2; 1 4])
@test_throws ErrorException compute_MM(my_test_NEP, [1 2; 1 4], diagm(0 => [1,2]))
@test_throws ErrorException size(my_test_NEP,1)

end

@bench @testset "Types (sumnep)" begin
    A0 = ones(3, 3)
    A1 = Matrix(1.0I, 3, 3)
    nep1 = DEP([A0,A1])
    B0 = reverse(Matrix(1.0I, 3, 3), dims = 1)
    B1 = Matrix(3.0I, 3, 3) + ones(3,3)
    nep2 = PEP([B0,B1])
    λ = 1+1im
    sumnep = SumNEP(nep1, nep2)
    M = compute_Mder(sumnep, λ)
    M1 = compute_Mder(nep1, λ)
    M2 = compute_Mder(nep2, λ)
    @test (M1+M2) ≈ M
end

# compute_Mlincomb for PEPs
@bench @testset "compute_Mlincomb PEP" begin
    nep=nep_gallery("pep0");
    n=size(nep,1);
    Random.seed!(0);
    λ = rand()+rand()*im; V=randn(n,3); a=rand(3);
    # test against another way to compute Mlincomb
    z=compute_Mlincomb(nep,λ,V,a)
    z2=compute_Mlincomb_from_MM(nep,λ,V,a)
    @test norm(z-z2)<sqrt(eps())*100
    # test that the function compute_Mlincomb does not overwrite the input
    λ = 0; V=randn(n,3); W=copy(V); a=rand(3);
    z=compute_Mlincomb(nep,λ,V,a);
    # test against another way to compute Mlincomb
    z2=compute_Mlincomb_from_MM(nep,λ,V,a)
    @test norm(z-z2)<sqrt(eps())*100
end


# compute_Mlincomb for DEPs
@bench @testset "compute_Mlincomb DEP" begin
    Random.seed!(0);
    n=100; A1=rand(n,n); A2=rand(n,n); A3=rand(n,n);
    tau1=0; tau2=1.3; tau3=.1;
    nep=DEP([A1,A2,A3],[tau1,tau2,tau3])
    n=size(nep,1);
    λ = rand()+rand()*im; V=randn(n,3); a=rand(3);
    # test against another way to compute Mlincomb
    z=compute_Mlincomb(nep,λ,V,a)
    z2=compute_Mlincomb_from_MM(nep,λ,V,a)
    @test norm(z-z2)<sqrt(eps())*100
    # test that the function compute_Mlincomb does not overwrite the input
    λ = 0; V=randn(n,3); W=copy(V); a=rand(3);
    z=compute_Mlincomb(nep,λ,V,a);
    # test against another way to compute Mlincomb
    z2=compute_Mlincomb_from_MM(nep,λ,V,a)
    @test norm(z-z2)<sqrt(eps())*100
end

# compute_Mlincomb for DerSPMF
@bench @testset "compute_Mlincomb DerSPMF" begin
    A0=[1 3; 4 5]; A1=[3 4; 5 6];
    id_op=S -> one(S) # Note: We use one(S) to be valid both for matrices and scalars
    exp_op=S -> exp(S)
    nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
    m=5; σ1=3.5; σ2=9;
    Dnep=DerSPMF(nep,σ1,m)
    DDnep=DerSPMF(Dnep,σ2,m)
    Random.seed!(0);
    V=randn(2,5);
    @test norm( compute_Mlincomb(nep,0,V)-compute_Mlincomb(Dnep,0,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ1,V)-compute_Mlincomb(Dnep,σ1,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ2,V)-compute_Mlincomb(Dnep,σ2,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,0,V)-compute_Mlincomb(DDnep,0,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ1,V)-compute_Mlincomb(DDnep,σ1,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ2,V)-compute_Mlincomb(DDnep,σ2,V) )<sqrt(eps())*100

    n=10;
    nep=nep_gallery("dep0",n)
    m=5; σ1=rand()+rand()*im; σ2=rand()+rand()*im;
    Dnep=DerSPMF(nep,σ1,m)
    DDnep=DerSPMF(Dnep,σ2,m)
    V=randn(n,5);
    @test norm( compute_Mlincomb(nep,0,V)-compute_Mlincomb(Dnep,0,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ1,V)-compute_Mlincomb(Dnep,σ1,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ2,V)-compute_Mlincomb(Dnep,σ2,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,0,V)-compute_Mlincomb(DDnep,0,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ1,V)-compute_Mlincomb(DDnep,σ1,V) )<sqrt(eps())*100
    @test norm( compute_Mlincomb(nep,σ2,V)-compute_Mlincomb(DDnep,σ2,V) )<sqrt(eps())*100
end
