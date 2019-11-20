using NonlinearEigenproblems, LinearAlgebra, Random
#import NonlinearEigenproblems.build_CORKPencil



# test on the following DEP
# M(λ)=-λI+A0+exp(-λ)vv'
Random.seed!(0);

n=2;
A0=rand(n,n)
A0=[1.0 3.0; -1.0 2.0];
v=reshape([-1.0 ; 1]/sqrt(2),n,1);

#rk=1;
#v=rand(n,1)
#u=randn(n,1)
dep=DEP([A0,v*v'],[0,1]);
cp_org=compute_CORKPencil(dep,IarCorkLinearization(d=4))


cplr=lowRankCompress(cp_org,1,1);


(AA_LR,BB_LR)=buildPencil(cplr);
evps_lr=eigen(AA_LR,BB_LR).values;


AA,BB=buildPencil(cp_org);
evps=eigen(AA,BB).values;
#

evps_lr_small=evps_lr[findall(abs.(evps_lr) .<5)]
for i = 1:length(evps_lr_small)
    @show minimum(abs.(evps_lr_small[i] .- evps))
end


A0=[1.0 3.0; -1.0 2.0]/10;
v=reshape([-1.0 ; 1]/sqrt(2),n,1);

Av=[-A0-v*v']
Bv=[-one(A0)-v*v']
BvLR=[v/2, -v/3, v/4, -v/5, v/6, -v/7,  v/8, -v/9]
AvLR=zero.(BvLR);
Z=v;

#
d=9;
M=diagm( 0 =>  ones(d) )[2:end,:]
N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]

cplr2=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z);

(AA,BB)=buildPencil(cplr2);
λ=eigen(AA,BB).values[end];
minimum(svdvals(A0-λ*I+v*v'*exp(-λ)))

#
#d=length(cp_org.Av);
#dtilde=2;
#cp_org.Bv[dtilde+1]
#
#Z=svd(cp_org.Bv[dtilde+1]).V[:,1:rk];
#
#Bvtilde=map(i-> cp_org.Bv[i]*Z,  (dtilde+1):d);
#Avtilde=map(i-> cp_org.Av[i]*Z,  (dtilde+1):d);
#
### Note info not in CORK paper:  (same for N)
###    M11: (dtilde-1) x (dtilde)
###    M21: (d-dtilde) x (dtilde)
###    M22: (d-dtilde) x (d-dtilde)
#M11=cp_org.M[1:(dtilde-1),1:dtilde];
#M21=cp_org.M[(dtilde):end,1:dtilde]
#M22=cp_org.M[(dtilde):end,(dtilde+1):end]
#
#N11=cp_org.N[1:(dtilde-1),1:dtilde];
#N21=cp_org.N[(dtilde):end,1:dtilde]
#N22=cp_org.N[(dtilde):end,(dtilde+1):end]
#
#In=Matrix{Float64}(I,n,n);
#Idtilde=Matrix{Float64}(I,rk,rk);
#
#Btilde1=[hcat(cp_org.Bv[1:dtilde]...) hcat(Bvtilde...) ]
#Btilde2=[kron(N11,In) zeros((dtilde-1)*n,(d-dtilde)*rk)]
#Btilde3=[kron(N21,Z') kron(N22,Idtilde)];
#Btilde=vcat(Btilde1,Btilde2,Btilde3);
#
#
#Atilde1=[hcat(cp_org.Av[1:dtilde]...) hcat(Avtilde...) ]
#Atilde2=[kron(M11,In) zeros((dtilde-1)*n,(d-dtilde)*rk)]
#Atilde3=[kron(M21,Z') kron(M22,Idtilde)];
#Atilde=vcat(Atilde1,Atilde2,Atilde3);
#
#
#
#
#evps1=eigen(Atilde,Btilde);
#AA,BB=build_CORKPencil(cp_org);
#evps2=eigen(AA,BB);
#
#minimum(abs.(evps1.values[1] .- evps2.values))



#zeros(rk*(d-dtilde-1),  )

#Avtilde=

#
#
#Av=[A0]
#Bv=[-one(A0)-v*v']
#BvLR=[v/2, v/3, v/4, v/5]
#AvLR=[zero(v),zero(v),zero(v),zero(v)]
#Z=v;
#
#d=4;
#M=diagm( 0 =>  ones(d) )[2:end,:]
#N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]
#
#c=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z)
#
## hardcoded contruction of the linearization
#dd=length(Av)
#M11=M[1:dd,1:dd]; 			N11=N[1:dd,1:dd];
#M21=M[dd+1:end,1:dd];		N21=N[dd+1:end,1:dd];
#M22=M[dd+1:end,dd+1:end];	N22=N[dd+1:end,dd+1:end];
#O=kron(M[1:dd,dd+1:end],zeros(n,d-dd));
#II1=one(A0);	II2=[1]; # II2=one(Z)
#AA=	[hcat(Av...) hcat(AvLR...);
#	[kron(M11,II1) 			zeros(size(kron(M11,II1),1),size(kron(M22,II2),2));
#	 kron(M21,copy(Z')) 	kron(M22,II2)]]
##BB=[hcat(Bv...) hcat(BvLR...);]
#
##AA=[hcat(Av...) hcat(AvLR...); kron(M,II)], [hcat(Bv...); kron(N,II)]
#
