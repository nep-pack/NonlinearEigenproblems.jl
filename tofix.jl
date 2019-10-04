using NonlinearEigenproblems, Random, SparseArrays, LinearAlgebra, Revise
import ..NEPSolver.ilan;
include("src/method_ilan.jl");


nep=nep_gallery("dep_symm_double",10);
v0=ones(size(nep,1));
λ,v=ilan(nep;v=v0,tol=1e-5,neigs=3,logger=1,check_error_every=1);
