using NonlinearEigenproblems, Random, SparseArrays, Revise, BenchmarkTools
import ..NEPSolver.ilan;
import ..NEPSolver.tiar;

include("../src/method_ilan.jl");
include("../src/method_tiar.jl");


# load the Voss symmetric DEP
n=50; nep=nep_gallery("dep_symm_double",n)
#@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=Inf) end
#@btime begin iar(Float64,nep;Neig=200,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=Inf) end

mm=400
v0=rand(Float64,n^2)
ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0)
@time begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0)
@time begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
1
