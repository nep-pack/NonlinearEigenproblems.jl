using NonlinearEigenproblems, Random, SparseArrays, Revise, BenchmarkTools
import ..NEPSolver.ilan;
import ..NEPSolver.tiar;

include("../src/method_ilan.jl");
include("../src/method_tiar.jl");


# load the Voss symmetric DEP
n=100; nep=nep_gallery("dep_symm_double",n)

# n=100000
# Random.seed!(1) # reset the random seed
# K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
# A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
# A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
# nep=DEP([A1,A2],[0,1])


mm=200
v0=rand(Float64,n^2)
@time begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@time begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=300
v0=rand(Float64,n^2)
@time begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@time begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=600
v0=rand(Float64,n^2)
@time begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@time begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

1
