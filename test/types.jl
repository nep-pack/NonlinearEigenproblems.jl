# Tests for core functionality

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems
using Test

@bench @testset "NEW" begin

n=3;
A1 = rand(n,n);
A2 = rand(n,n);
A3 = rand(n,n);
A4 = rand(n,n);

f1 = S -> one(S)
f2 = S -> -S
f3 = S -> exp(-S)
f4 = S -> sqrt(S)
nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4]);

## Mlincomb_tests
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


end
