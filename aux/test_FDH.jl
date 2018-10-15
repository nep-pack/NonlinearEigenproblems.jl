using NonlinearEigenproblems, Random, SparseArrays, LinearAlgebra

n=100;
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2))
A2 = sparse(K, J, rand(3*n-2))
A3 = sparse(K, J, rand(3*n-2))
A1 = A1+A1';
A2 = A2+A2';
A3 = A3+A3';


f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)

nep=SPMF_NEP([A1,A2,A3],[f1,f2,f3])

fv=get_fv(nep); p=length(fv)
Av=get_Av(nep)

T=ComplexF64
m=10
μ=0
S=diagm(0 => μ*ones(T,2m)) + diagm(-1 => (1:2m-1))
# derivative matrix
fD=zeros(2*m,p)
for t=1:p fD[:,t]=fv[t](S)[:,1] end
FDH=Vector{Array{T,2}}(undef,p)
for t=1:p FDH[t]=zeros(T,m,m)
    for i=1:m for j=1:m
        FDH[t][i,j]=fD[i+j,t];
    end end
end
