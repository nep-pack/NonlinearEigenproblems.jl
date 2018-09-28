# Tests for SPMF-code

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random
using SparseArrays

@testset "SPMF stability" begin
    n=5;
    Random.seed!(1)
    A0=sparse(randn(n,n));
    A1=sparse(randn(n,n));
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
