# CORK data types



using NonlinearEigenproblems
using Test,Random
using LinearAlgebra


@bench @testset "CORK types" begin

    # Construct a low-rank NEP
    A0=[1.0 3.0; -1.0 2.0];
    v=reshape([-1.0 ; 1]/sqrt(2),2,1);
    dep=DEP([A0,v*v'],[0,1]);

    # Form a CORKPencil
    cp_org=CORKPencil(dep,IarCorkLinearization(d=9))
    # .. and compress it
    cplr=lowRankCompress(cp_org,1,1);
    (AA_LR,BB_LR)=buildPencil(cplr);

    # Manually construct a pencil
    Av=[-A0-v*v']
    Bv=[-one(A0)-v*v']
    BvLR=[v/2, -v/3, v/4, -v/5, v/6, -v/7,  v/8, -v/9]
    AvLR=zero.(BvLR);
    Z=v;

    d=9;
    M=diagm( 0 =>  ones(d) )[2:end,:]
    N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]
    cplr2=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z);

    (AA2,BB2)=buildPencil(cplr2);

    # ... and check that it's the same as we had before
    @test AA2 ≈ AA_LR
    @test BB2 ≈ BB_LR

end
