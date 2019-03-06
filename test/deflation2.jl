# Run tests for the deflation

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra
using BlockArrays


function create_spmf_dnep(nep::AbstractSPMF,S0,V0)
    T=eltype(V0);
    Av_org=get_Av(nep);
    fv_org=get_fv(nep);
    m=size(fv_org,1);
    p=size(V0,2);
    n0=size(nep,1);

    m1=m;
    m2=m*p+1;


    # spmf1: Create the "old" part
    A1=Vector{Matrix}(undef,m1);
    for k=1:m
        A1[k]=zeros(T,n+p,n+p)
        A1[k][1:n,1:n]=Av_org[k];
    end
    spmf1=SPMF_NEP(A1,fv_org,check_consistency=false)
    # spmf2: Create the additional deflation terms:
    A2=Vector{BlockMatrix{T}}(undef,m2);
    L2=Vector{Matrix{T}}(undef,m2);
    U2=Vector{Matrix{T}}(undef,m2);
    fv2=Vector{Function}(undef,m2);
    (λ,X)=eigen(S0);
    count=0;
    for i=1:p
        ei=zeros(p); ei[i]=1;
        @show ei
        y=(V0*(X*ei));
        x=(ei'/X);
        for r=1:m
            count=count+1;
            @show m+count
            A2[count]=BlockMatrix{T}(undef_blocks,[n0,p],[n0,p]);
            A2[count].blocks[1,1]=zeros(n,n);
            A2[count].blocks[1,2]=(Av_org[r]*y)*x;
            A2[count].blocks[2,1]=zeros(p,n);
            A2[count].blocks[2,2]=zeros(p,p);
            L2[count] = reshape([(Av_org[r]*y) ;zeros(p)],n+p,1);
            U2[count]=reshape([zeros(n);x'],n+p,1);
            fv2[count]=S-> inv(S-λ[i]*one(S))*fv_org[r](S);
        end
    end
    A2[end]=BlockMatrix{T}(undef_blocks,[n0,p],[n0,p]);
    A2[end].blocks[1,1]=zeros(n,n);
    A2[end].blocks[1,2]=zeros(n,p);
    A2[end].blocks[2,1]=Matrix(V0');
    A2[end].blocks[2,2]=zeros(p,p);
    L2[end]=[zeros(n,p);Matrix{T}(I,p,p)]
    U2[end]=[Matrix(V0);zeros(p,p)]
    fv2[end]= S->one(S);
    spmf2=LowRankFactorizedNEP(L2,U2,fv2);
    spmf2b=SPMF_NEP(A2,fv2,check_consistency=false);
    z=norm(compute_Mder(spmf2,4.1+1im)-compute_Mder(spmf2b,4.1+1im))
    @show z

    return SumNEP(spmf1,spmf2);
end


nep=nep_gallery("dep0");
n=size(nep,1);
(λ,v)=newton(nep,v=ones(n),λ=0.8)
v=v/norm(v);
S0=reshape([λ],1,1);
V0=reshape(v,n,1);



dnep=effenberger_deflation(nep,S0,V0)

dnep2=create_spmf_dnep(nep,S0,V0)


W=randn(size(dnep,1),3);
#W[3:end].=0;
Z=randn(3,3);
@show norm(compute_MM(dnep,Z,W)-compute_MM(dnep2,Z,W))

(λ2,v2)=newton(dnep2,λ=-0.3,armijo_factor=0.9,displaylevel=1,maxit=100,v=ones(n+1))
v2=v2/norm(v2);



#Z=randn(1,1); X=randn(n+1,1);
#compute_MM(dnep,Z,X)-compute_MM(dnep2,Z,X)

# Create the new invariant pair
V=V0;
S=S0;

V1=zeros(n,size(V,2)+1);
S1=zeros(size(S,1)+1,size(S,2)+1);
V1[1:n,1:end-1]=V[1:n,:];
V1[1:n,end]=v2[1:n];
S1[1:end-1,1:end-1]=S;
S1[1:end,end]=[v2[n+1:end];λ2];


dnep3=create_spmf_dnep(nep,S1,V1)
dnep3b=effenberger_deflation(nep,S1,V1);
Z=randn(1,1); X=randn(n+2,1);
compute_MM(dnep3b,Z,X)-compute_MM(dnep3,Z,X)

(λ3,v3)=newton(dnep3,λ=1.4im,armijo_factor=0.9,displaylevel=1,maxit=100,v=ones(n+2))
v3=v3/norm(v3);


V=V1;
S=S1;

V1=zeros(ComplexF64,n,size(V,2)+1);
S1=zeros(ComplexF64,size(S,1)+1,size(S,2)+1);
V1[1:n,1:end-1]=V[1:n,:];
V1[1:n,end]=v3[1:n];
S1[1:end-1,1:end-1]=S;
S1[1:end,end]=[v3[(n+1):end];λ3];

compute_MM(nep,S1,V1)
