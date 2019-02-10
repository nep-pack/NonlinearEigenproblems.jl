# Test for nep helper creators
using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test

@bench @testset "NEPTypes: creat help" begin
    function my_Mder(s,der)
        A0=ones(Float64,3,3);
        A1=ones(Float64,3,3)*3; A1[2,3]=0;
        return (3^der)*exp(3*s)*A0+factorial(der)*((-1)^der)*A1/(3.3+s)^(der+1)
    end
    AA0=ones(Float64,3,3);
    AA1=ones(Float64,3,3)*3; AA1[2,3]=0;
    nep_ref=SPMF_NEP([AA0,AA1],[s->exp(3*s), S->inv(3.3*I+S)]);
    nep=Mder_NEP(my_Mder,3)

    λ=3.3;
    norm(compute_Mder(nep,λ,4)-compute_Mder(nep_ref,λ,4))

    λ=-1.2+0.2im;
    X=[ -0.715845   0.865534   0.254796 ; -0.856405  -0.482516   0.0265129;
  1.02593   -0.62892   -2.09615  ]
    norm(compute_Mlincomb(nep,λ,X)-compute_Mlincomb(nep_ref,λ,X))
end
