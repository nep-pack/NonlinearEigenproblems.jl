using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra
import ..NEPTypes.AbstractSPMF

struct DerSPMF
    spmf::SPMF_NEP
    fD::Matrix
    σ::Number
end

## one constructor takes spmf as input and compute the derivatives
function DerSPMF(spmf::SPMF_NEP,σ::Number)
      # Compute DD-matrix from get_fv(spmf)
      m=100;    # it should be an input
      p=length(nep.fi)
      fD=zeros(m,p)
      # matrix for the computation of derivatives
      SS=diagm(0 => σ*ones(2m+2)) + diagm(-1 => (1:2m+1))
      fD=zeros(2*m+2,p)
      for t=1:p fD[:,t]=nep.fi[t](SS)[:,1] end
      return DerSPMF(spmf,fD,σ);
end


n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
A4 = sparse(K, J, rand(3*n-2)); A4 = A4+A4';

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
f4= S -> exp(-S)

nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
σ=rand()
DD=rand(2,2)
Dnep=DerSPMF(nep,DD,σ)
Dnep2=DerSPMF(nep,σ)
