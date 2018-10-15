using NonlinearEigenproblems, Random, SparseArrays, Revise
import ..NEPTypes.AbstractSPMF

struct DerSPMF
    spmf::SPMF_NEP
    DD::Matrix
    σ::Number
end

## one constructor takes spmf as input and compute the derivatives
function DerSPMF(spmf::SPMF_NEP,σ::Number)
      # Compute Q-matrix from get_fv(spmf)
      DD=rand(2,2)
      return DerSPMF(spmf,DD,σ);
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
