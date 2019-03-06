# Run tests for the deflation

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra
using BlockArrays


function create_spmf_dnep(nep::AbstractSPMF,V0,S0)
    T=eltype(V0);
    Av_org=get_Av(nep);
    fv_org=get_fv(nep);
    m=size(fv_org,1);
    p=size(V0,2);
    n0=size(nep,1);
    m_new=m*(p+1)+1;
    A_new=Vector{BlockMatrix{T}}(undef,m_new);
    fv_new=Vector{Function}(undef,m_new);

    # Create the "old" part
    for k=1:m
        A_new[k]=BlockMatrix{T}(undef_blocks,[n0,p],[n0,p]);
        A_new[k].blocks[1,1]=Av_org[k];
        A_new[k].blocks[1,2]=zeros(n,p);
        A_new[k].blocks[2,1]=zeros(p,n);
        A_new[k].blocks[2,2]=zeros(p,p);
        fv_new[k]=fv_org[k];
    end

    # Create the additional deflation terms:
    (λ,X)=eigen(S0);
    count=0;
    @show size(A_new)
    for i=1:p
        ei=zeros(p); ei[i]=1;
        @show ei
        y=(V0*(X\ei));
        x=(ei'*X);
        for r=1:m
            count=count+1;
            @show m+count
            A_new[m+count]=BlockMatrix{T}(undef_blocks,[n0,p],[n0,p]);
            A_new[m+count].blocks[1,1]=zeros(n,n);
            A_new[m+count].blocks[1,2]=(Av_org[r]*y)*x;
            A_new[m+count].blocks[2,1]=zeros(p,n);
            A_new[m+count].blocks[2,2]=zeros(p,p);
            fv_new[m+count]=S-> -inv(λ[i]*I-S)*fv_org[r](S);
        end
    end
    A_new[end]=BlockMatrix{T}(undef_blocks,[n0,p],[n0,p]);
    A_new[end].blocks[1,1]=zeros(n,n);
    A_new[end].blocks[1,2]=zeros(n,p);
    A_new[end].blocks[2,1]=Matrix(V0');
    A_new[end].blocks[2,2]=zeros(p,p);
    fv_new[end]= S->one(S);
    return SPMF_NEP(A_new,fv_new,check_consistency=false);
end


nep=nep_gallery("dep0");
n=size(nep,1);
(λ,v)=newton(nep,v=ones(n))
v=v/norm(v);
S0=reshape([λ],1,1);
V0=reshape(v,n,1);



dnep=effenberger_deflation(nep,S0,V0)

dnep2=create_spmf_dnep(nep,V0,S0)


W=randn(size(dnep,1),3);
#W[3:end].=0;
Z=randn(3,3);
norm(compute_MM(dnep,Z,W)-compute_MM(dnep2,Z,W))

(λ2,v2)=newton(dnep2,λ=3,armijo_factor=0.9,displaylevel=1,maxit=100)
