using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra, BenchmarkTools

n=1000000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
A4 = sparse(K, J, rand(3*n-2)); A4 = A4+A4';

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
f4= S -> sqrt(10*one(S)-S)

nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
σ=rand()
DD=rand(2,2)
#Dnep=DerSPMF(nep,DD,σ)
Dnep=DerSPMF(nep,σ,100)
m=40
V=rand(n,m)




function compare()
      @btime begin z1=compute_Mlincomb(nep,σ,V) end
      @btime begin z2=compute_Mlincomb(Dnep,σ,V) end
end

compare()
z1=compute_Mlincomb(nep,σ,V)
z2=compute_Mlincomb(Dnep,σ,V)
norm(z1-z2)

σ2=rand()
Dnep2=DerSPMF(Dnep,σ2,100)

function compare2()
      @btime begin z1=compute_Mlincomb(nep,σ2,V) end
      @btime begin z2=compute_Mlincomb(Dnep,σ2,V) end
      @btime begin z3=compute_Mlincomb(Dnep2,σ2,V) end
end
compare2()

z1=compute_Mlincomb(nep,σ2,V)
z2=compute_Mlincomb(Dnep,σ2,V)
z3=compute_Mlincomb(Dnep2,σ2,V)
display(norm(z1-z2))
display(norm(z1-z3))
