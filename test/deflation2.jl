# Run tests for the deflation

using NonlinearEigenproblems
using Test
using LinearAlgebra
using SparseArrays;
using BenchmarkTools;




#
#@show "nlevp_native"
#nep=nep_gallery("nlevp_native_gun");
#n=size(nep,1);
#(λ,v)=augnewton(nep,v=ones(n),λ=200^2,tol=1e-11)
#
#dnep1=deflate_eigpair(nep,λ,v,mode=:Generic);
#dnep2=deflate_eigpair(nep,λ,v,mode=:SPMF);
#dnep3=deflate_eigpair(nep,λ,v,mode=:MM);
#
#@btime augnewton(dnep1,v=ones(size(dnep1,1)),λ=250^2,
#                      tol=1e-11,maxit=300,armijo_factor=0.5)
#@btime augnewton(dnep2,v=ones(size(dnep1,1)),λ=250^2,
#                      tol=1e-11,maxit=300,armijo_factor=0.5)
#@btime augnewton(dnep3,v=ones(size(dnep1,1)),λ=complex(250^2),
#                      tol=1e-11,maxit=300,armijo_factor=0.5)
#asd()

nep=nep_gallery("dep0_sparse",100);
n=size(nep,1);
local λ,v;
(λ,v)=quasinewton(nep,v=ones(n),λ=2,tol=1e-11,logger=1,maxit=200)
dnep1=deflate_eigpair(nep,λ,v,mode=:Generic);
@btime mslp(dnep1,λ=1im,tol=1e-11)
dnep2=deflate_eigpair(nep,λ,v,mode=:SPMF);
@btime mslp(dnep2,λ=1im,tol=1e-11)
dnep3=deflate_eigpair(nep,λ,v,mode=:MM);
@btime mslp(dnep3,λ=1im,tol=1e-11)


asd()
begin
  nep=nep_gallery("dep0");
  n=size(nep,1);
  local λ,v;
  (λ,v)=quasinewton(nep,v=ones(n),λ=0,tol=1e-11)
  local dnep
  dnep=deflate_eigpair(nep,λ,v,mode=:Generic);
  for k=1:4
      (λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=-1+0.1im,
                      tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)

      (λ2,V2)=get_deflated_eigpairs(dnep,λ,v)
      @show norm(compute_Mlincomb(nep,λ2[end],V2[:,end]))
      dnep=deflate_eigpair(dnep,λ,v);
  end
  (λ,V)=get_deflated_eigpairs(dnep)
  @show norm(compute_Mlincomb(nep,λ[end],V[:,end]))
  @show norm(compute_Mlincomb(nep,λ[1],V[:,1]))
end

begin
  nep=nep_gallery("nlevp_native_gun");
  n=size(nep,1);
  local λ,v;
  (λ,v)=augnewton(nep,v=ones(n),λ=200^2,tol=1e-11)
  local dnep
  dnep=deflate_eigpair(nep,λ,v,mode=:Generic);
  for k=1:2
      (λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=200^2,
                      tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)

      (λ2,V2)=get_deflated_eigpairs(dnep,λ,v)
      @show norm(compute_Mlincomb(nep,λ2[end],V2[:,end]))
      dnep=deflate_eigpair(dnep,λ,v);
  end
  (λ,V)=get_deflated_eigpairs(dnep)
  @show norm(compute_Mlincomb(nep,λ[end],V[:,end]))
  @show norm(compute_Mlincomb(nep,λ[1],V[:,1]))
end


asd()


begin
  nep=nep_gallery("dep0");
  n=size(nep,1);
  local λ,v;
  (λ,v)=quasinewton(nep,v=ones(n),λ=0,tol=1e-11)
  local dnep
  dnep=deflate_eigpair(nep,λ,v,mode=:MM);
  for k=1:4
      (λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=-1+0.1im,
                      tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)

      (λ2,V2)=get_deflated_eigpairs(dnep,λ,v)
      @show norm(compute_Mlincomb(nep,λ2[end],V2[:,end]))
      dnep=deflate_eigpair(dnep,λ,v);
  end
  (λ,V)=get_deflated_eigpairs(dnep)
  @show norm(compute_Mlincomb(nep,λ[end],V[:,end]))
  @show norm(compute_Mlincomb(nep,λ[1],V[:,1]))
end

#(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=0,
#                tol=1e-11,displaylevel=1,maxit=300,armijo_factor=0.5)
#dnep=deflate_eigpair(dnep,λ,v);
#(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=0.1im	,
#                tol=1e-11,displaylevel=1,maxit=300,armijo_factor=0.5)
#dnep=deflate_eigpair(dnep,λ,v);
#(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=0.1im	,
#                tol=1e-11,displaylevel=1,maxit=300,armijo_factor=0.5)
#dnep=deflate_eigpair(dnep,λ,v);

asd()

#nep=nep_gallery("nlevp_native_gun");
n=size(nep,1);
(λ,v)=quasinewton(nep,v=ones(n),λ=150^2,tol=1e-11)
dnep=deflate_eigpair(nep,λ,v,mode=:Generic);

asd()
(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)

dnep=deflate_eigpair(dnep,λ,v);
(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)
dnep=deflate_eigpair(dnep,λ,v);
(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)

(λ,V)=get_deflated_eigpairs(dnep)

asd()

nep=nep_gallery("nlevp_native_gun");
n=size(nep,1);
(λ,v)=quasinewton(nep,v=ones(n),λ=150^2,tol=1e-11)

dnep=deflate(dnep,λ,v)

(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)
dnep=deflate(dnep,λ,v)

(λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)



asd()

# Wanted usage:
# (λ,v)=quasinewton(nep)
# dnep=deflate(nep,λ,v)
# # dnep=deflate(nep,λ,v,mode=:SPMF); # Default
# # dnep=deflate(nep,λ,v,mode=:NaiveSPMF)
# # dnep=deflate(nep,λ,v,mode=:Generic)
# (λ,v)=quasinewton(dnep)
# dnep=deflate(dnep,λ,v)
# # dnep=deflate_schur(nep,S,V)
#
# # Extraction usage 1:
# (D,V)=deflated_eigpairs(dnep)
# # Extraction usage 2:
# (D,V)=deflated_eigpairs(dnep,λ,v)

# # Extraction Schur pairs
# (S,X)=deflated_schurfact(dnep)
# (S,X)=deflated_schurfact(dnep,λ,v)
#

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

    @show 1
    for k=1:3
        normalize_schur_pair!(S1,V1);
        @show norm(compute_MM(nep,S1,V1))
        dnep=create_spmf_dnep(nep,S1,V1)
        #λ,V=get_deflated_eigpairs(dnep);
        @show norm(compute_Mlincomb(nep,λ[end],V[:,end]))
        (λ,v)=augnewton(dnep,v=ones(size(dnep,1)),λ=300^2,
                        tol=1e-11,logger=1,maxit=300,armijo_factor=0.5)

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

(λ2,v2)=resinv(dnep1_new,λ=290^2,armijo_factor=0.9,logger=1,maxit=100,v=ones(n+1),tol=1e-12)
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
                    logger=1,maxit=100,v=ones(n+2),tol=1e-12)

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
                    logger=1,maxit=200,v=ones(n+3),tol=1e-12)

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
                    logger=1,maxit=200,v=ones(n+4),tol=1e-12)

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
