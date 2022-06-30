# Tests for SPMF-code

using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random
using SparseArrays

@testset "SPMF stability" begin
    n=5;
    Random.seed!(1)
    A0=sparse([  0.297288    0.311111   0.583708  -0.0504512  -1.72976
  0.382396    2.29509    0.963272  -0.693654    0.795949
 -0.597634   -2.26709    0.458791  -1.77383     0.670062
 -0.0104452   0.529966  -0.522337   0.120766    0.550852
          -0.839027    0.431422   0.408396  -0.757633   -0.0633746]);
    A1=sparse([  1.33694    -0.165136   0.763499  -1.28905     1.29992
 -0.0731486  -2.11537    0.378599  -0.33527     0.206364
 -0.745464   -0.066768  -0.645597   0.0704676  -1.00886
 -1.22006     1.22246   -0.664646   0.341794   -0.850056
 -0.0531773   0.567695  -1.80303    1.73517     1.12941 ]);

    t=3.0

    minusop = S -> -S
    oneop = S -> one(S)
    expmop = S -> exp(-t*S)
    fi=[minusop, oneop, expmop];

    J = SparseMatrixCSC(1.0I, n, n)
    nep1 = SPMF_NEP([J, A0, A1], fi)
    nep2 = SPMF_NEP([J, A0, A1], fi; Schur_fact=true)

    m=4;
    V=randn(n,m);
    λ=3;
    S=randn(m,m);
    @inferred compute_MM(nep1,S,V)

    compute_Mlincomb(nep1,λ,V)


    @inferred compute_Mlincomb(nep1,3,V)
    @inferred compute_Mlincomb(nep1,3+1im,V)

    nep2=SPMF_NEP(get_Av(nep1),get_fv(nep1),Ftype=Float64);

    @inferred compute_Mlincomb(nep1,3,V)
    @inferred compute_Mlincomb(nep1,3+1im,V)

end
