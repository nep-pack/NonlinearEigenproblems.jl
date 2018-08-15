# Tests for core functionality

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using IterativeSolvers
    using Base.Test
end


struct TestNEP <: NEP
       dummy
end

@testset "Core" begin

nep=nep_gallery("dep0");
n=size(nep,1);

## Mlincomb_tests
V=complex(randn(n,3));
λ=0.3+1im;
T=typeof(λ);
z1=compute_Mlincomb(nep,λ,V)
z2=compute_Mlincomb(nep,λ,V,a=[T(1),T(1),T(1)])
@test z1==z2
z3=compute_Mder(nep,λ,0)*V[:,1]+ compute_Mder(nep,λ,1)*V[:,2]+ compute_Mder(nep,λ,2)*V[:,3]
@test z1 ≈ z3

## Test Mlincomb starting counting with kth derivative
z2=compute_Mlincomb(nep,λ,V,a=[T(0),T(0),T(1)])
z3=compute_Mder(nep,λ,2)*V[:,3]
@test z2≈z3
z4=compute_Mlincomb(nep,λ,V[:,3],[T(1)], 2);
@test z3≈z4
#

## Simple coverage. Test that methods throw error if not defined
my_test_NEP = TestNEP([1 2; 3 4])
@test_throws ErrorException compute_Mder(my_test_NEP, 1+1im, 2)
@test_throws ErrorException compute_Mlincomb(my_test_NEP , 1+1im, [1 2; 1 4])
@test_throws ErrorException compute_MM(my_test_NEP,[1 2; 1 4],diagm([1,2]))
@test_throws ErrorException size(my_test_NEP,1)

end


@testset "Types (sumnep)" begin
    A0=ones(3,3);
    A1=eye(3,3);
    nep1=DEP([A0,A1])
    B0=flipdim(eye(3,3),1)
    B1=3*eye(3,3)+ones(3,3);
    nep2=PEP([B0,B1]);
    λ=1+1im;
    sumnep=SumNEP(nep1,nep2);
    M=compute_Mder(sumnep,λ);
    M1=compute_Mder(nep1,λ);
    M2=compute_Mder(nep2,λ);
    @test (M1+M2)≈M
end
