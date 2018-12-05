using LinearAlgebra, Revise, GenericSVD, DelimitedFiles
m=200
# precompute the symmetrizer coefficients
setprecision(BigFloat,500)
G=zeros(BigFloat,m+1,m+1);
for i=1:m+1 G[i,1]=1/i end
for j=1:m
    for i=1:m+1
        G[i,j+1]=(G[i,j]*j)/(i+j);
    end
end
tolG=1e-14
U,S,V=svd(G);
precision(BigFloat)

SS=zeros(BigFloat,m+1,2)
SS[1:m+1,1]=1:m+1
SS[1:m+1,2]=S
SS=SS[1:50,:]
writedlm("G_svd200.csv",SS,",")
