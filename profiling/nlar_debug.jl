using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra, BenchmarkTools

nep=nep_gallery("dep0",100);
v0=ones(size(nep,1));
try
    λ,v=iar(nep;v=v0,tol=1e-5,neigs=10,maxit=10);
catch e
    λ=e.λ
    v=e.v
end
