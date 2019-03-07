# Run tests for the deflation

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra
using BlockArrays
using SparseArrays;

function create_spmf_dnep(nep::AbstractSPMF,S0,V0)
    Av_org=get_Av(nep);
    fv_org=get_fv(nep);
    m=size(fv_org,1);
    p=size(V0,2);
    n0=size(nep,1);



    m1=m;     # size of "old" part
    m2=m*p+1; # size of "deflation" part

    # spmf1: Create the "old" part
    A1=Vector{eltype(Av_org)}(undef,m1);
    for k=1:m
        A0k=Av_org[k];
        if (eltype(A1) <:  SparseMatrixCSC)
            (II,JJ,VV)=findnz(A0k)
            A1[k]=sparse(II,JJ,VV,n+p,n+p);
        else
            A1[k]=zeros(eltype(A0k),n+p,n+p)
            A1[k][1:n,1:n]=A0k;
        end
    end
    spmf1=SPMF_NEP(A1,fv_org,check_consistency=false)
    # spmf2: Create the additional deflation terms:

    #    Promote rules for the eltype:
    #    We may need to increase the eltype type size:
    T=promote_type(eltype(V0),eltype(S0),eltype(Av_org[1]));
    local T_LowRankFactor;
    if (eltype(Av_org) <: SparseMatrixCSC)
        T_LowRankFactor=SparseMatrixCSC{T,Int64};
    else
        T_LowRankFactor=Matrix{T};
    end
    L2=Vector{T_LowRankFactor}(undef,m2);
    U2=Vector{T_LowRankFactor}(undef,m2);
    fv2=Vector{Function}(undef,m2);
    (λtmp,X)=eigen(S0);
    λ::Vector{T}=λtmp[:]; # Ensure type
    count=0;
    for i=1:p
        ei=zeros(p); ei[i]=1;
        y=(V0*(X[:,i]));
        x=(ei'/X);
        for r=1:m
            count=count+1;
            # This will automatically convert to sparse / full
            L2[count] = reshape([(Av_org[r]*y) ;zeros(p)],n+p,1);
            U2[count] = reshape([zeros(n);x'],n+p,1);
            fv2[count]=S-> (S-λ[i]*one(S))\fv_org[r](S);
        end
    end
    # The constant term
    L2[m*p+1]=[zeros(n,p);Matrix{T}(I,p,p)]
    U2[m*p+1]=[Matrix(V0);zeros(p,p)]
    fv2[m*p+1]= S->one(S);
    spmf2=LowRankFactorizedNEP(L2,U2,fv2);

    return SumNEP(spmf1,spmf2);
end

function normalize_schur_pair!(S,V)
    (QQ,RR)=qr(V);
    V[:]=Matrix(QQ); # Use skinny QR-factorization
    S[:]=(RR*S)/RR;
end



# Wanted usage:
# (λ,v)=quasinewton(nep)
# dnep=deflate(nep,λ,v)
# # dnep=deflate(nep,λ,v,mode=:SPMF)
# # dnep=deflate(nep,λ,v,mode=:NaiveSPMF)
# # dnep=deflate(nep,λ,v,mode=:Generic)
# (λ,v)=quasinewton(dnep)
# dnep=deflate(dnep,λ,v)
# # dnep=deflate_schur(nep,S,V)
#
# # Extraction usage 1:
# (λ,v)=deflated_eigpairs(dnep,λ,v)
#
# # Extraction usage 2:
# (λ,v)=deflated_eigpairs(dnep)

nep=nep_gallery("nlevp_native_gun");
n=size(nep,1);
(λ,v)=quasinewton(nep,v=ones(n),λ=150^2,tol=1e-11)

#nep=nep_gallery("dep0_sparse");
#n=size(nep,1);
#(λ,v)=augnewton(nep,v=ones(n),λ=-0.4+0.3im,tol=1e-11)
v=v/norm(v);
S0=reshape([λ],1,1);
V0=reshape(v,n,1);
begin
    local S1=reshape([λ],1,1);
    local V1=reshape(v,n,1);

    for k=1:8
        normalize_schur_pair!(S1,V1);
        @show norm(compute_MM(nep,S1,V1))
        dnep=create_spmf_dnep(nep,S1,V1)
        (λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                        tol=1e-11,displaylevel=1,maxit=200,armijo_factor=0.5)

        V=V1;
        S=S1;
        V1=zeros(eltype(v),n,size(V,2)+1);
        S1=zeros(eltype(v),size(S,1)+1,size(S,2)+1);
        V1[1:n,1:end-1]=V[1:n,:];
        V1[1:n,end]=v[1:n];
        S1[1:end-1,1:end-1]=S;
        S1[1:end,end]=[v[n+1:end];λ];

        normalize_schur_pair!(S1,V1);

    end
    @show diag(S1)
    @show diag(V1'*V1)
end

asd()

# dnep1_eff=effenberger_deflation(nep,S0,V0)

dnep1_new=create_spmf_dnep(nep,S0,V0)

#W=randn(size(nep,1)+1,3);
##W[3:end].=0;
#Z=randn(3,3);
#@show norm(compute_MM(dnep,Z,W)-compute_MM(dnep2,Z,W))

(λ2,v2)=resinv(dnep1_new,λ=290^2,armijo_factor=0.9,displaylevel=1,maxit=100,v=ones(n+1),tol=1e-12)
v2=v2/norm(v2);


#Z=randn(1,1); X=randn(n+1,1);
#compute_MM(dnep,Z,X)-compute_MM(dnep2,Z,X)

# Create the new invariant pair
V=V0;
S=S0;

V1=zeros(eltype(v2),n,size(V,2)+1);
S1=zeros(eltype(v2),size(S,1)+1,size(S,2)+1);
V1[1:n,1:end-1]=V[1:n,:];
V1[1:n,end]=v2[1:n];
S1[1:end-1,1:end-1]=S;
S1[1:end,end]=[v2[n+1:end];λ2];

normalize_schur_pair!(S1,V1);

@show norm(compute_MM(nep,S1,V1))





dnep2_new=create_spmf_dnep(nep,S1,V1)
#dnep2_eff=effenberger_deflation(nep,S1,V1);
#Z=randn(1,1); X=randn(n+2,1);
#compute_MM(dnep3b,Z,X)-compute_MM(dnep3,Z,X)

(λ3,v3)=quasinewton(dnep2_new,λ=300^2,armijo_factor=0.9,
                    displaylevel=1,maxit=100,v=ones(n+2),tol=1e-12)

V=V1;
S=S1;

V1=zeros(ComplexF64,n,size(V,2)+1);
S1=zeros(ComplexF64,size(S,1)+1,size(S,2)+1);
V1[1:n,1:end-1]=V[1:n,:];
V1[1:n,end]=v3[1:n];
S1[1:end-1,1:end-1]=S;
S1[1:end,end]=[v3[(n+1):end];λ3];



@show norm(compute_MM(nep,S1,V1))


dnep3_new=create_spmf_dnep(nep,S1,V1)
#dnep2_eff=effenberger_deflation(nep,S1,V1);
#Z=randn(1,1); X=randn(n+2,1);
#compute_MM(dnep3b,Z,X)-compute_MM(dnep3,Z,X)

(λ4,v4)=quasinewton(dnep3_new,λ=295^2,armijo_factor=0.9,
                    displaylevel=1,maxit=200,v=ones(n+3),tol=1e-12)

normalize!(v4)

V=V1;
S=S1;

V1=zeros(ComplexF64,n,size(V,2)+1);
S1=zeros(ComplexF64,size(S,1)+1,size(S,2)+1);
V1[1:n,1:end-1]=V[1:n,:];
V1[1:n,end]=v4[1:n];
S1[1:end-1,1:end-1]=S;
S1[1:end,end]=[v4[(n+1):end];λ4];

normalize_schur_pair!(S1,V1);


@show norm(compute_MM(nep,S1,V1))




dnep4_new=create_spmf_dnep(nep,S1,V1)

(λ5,v5)=quasinewton(dnep4_new,λ=295^2,armijo_factor=0.9,
                    displaylevel=1,maxit=200,v=ones(n+4),tol=1e-12)

V=V1;
S=S1;

V1=zeros(ComplexF64,n,size(V,2)+1);
S1=zeros(ComplexF64,size(S,1)+1,size(S,2)+1);
V1[1:n,1:end-1]=V[1:n,:];
V1[1:n,end]=v5[1:n];
S1[1:end-1,1:end-1]=S;
S1[1:end,end]=[v5[(n+1):end];λ5];

normalize_schur_pair!(S1,V1);
#
#(QQ,RR)=qr(V1);
#V1=Matrix(QQ);
#S1=(RR*S1)/RR;
#

@show norm(compute_MM(nep,S1,V1))
