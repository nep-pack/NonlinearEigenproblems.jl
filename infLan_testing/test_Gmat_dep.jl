using LinearAlgebra, Revise
m=10
# precompute the symmetrizer coefficients
G=zeros(m+1,m+1);
for i=1:m+1 G[i,1]=1/i end
for j=1:m
    for i=1:m+1
        G[i,j+1]=(G[i,j]*j)/(i+j);
    end
end
tolG=1e-14
U,S,V=svd(G);
q=sum(S.>tolG*ones(length(S)))
U=U[:,1:q]; S=S[1:q]; V=V[:,1:q]
norm(G-U*diagm(0=>S)*V')


# precompute derivatives and FDH matrices
τ=2
f1 = S -> exp(-τ*S)
fv = [f1]

σ=2; γ=0.5
SS=diagm(0 => σ*ones(2m+2)) + diagm(-1 => γ*(1:2m+1))
p=1
fD=zeros(2*m+2,p)
for t=1:p fD[:,t]=fv[t](SS)[:,1] end

FDH=Vector{Array{Float64,2}}(undef,p)
for t=1:p
    FDH[t]=zeros(m+1,m+1)
    for i=1:m+1 for j=1:m+1
        FDH[t][i,j]=fD[i+j,t];
    end end
end
v=sqrt(τ*γ)*exp(-σ)*(-τ*γ).^(0:m)

norm(FDH[1]+v*v')
