using NonlinearEigenproblems, Random, SparseArrays, Revise, BenchmarkTools
import ..NEPSolver.ilan;
import ..NEPSolver.tiar;

include("../src/method_ilan.jl");
include("../src/method_tiar.jl");

# load the Voss symmetric DEP
n=20; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end




# load the Voss symmetric DEP
n=40; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end




# load the Voss symmetric DEP
n=80; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end



# load the Voss symmetric DEP
n=160; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end


# load the Voss symmetric DEP
n=320; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end


# load the Voss symmetric DEP
n=640; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end


# load the Voss symmetric DEP
n=1280; nep=nep_gallery("dep_symm_double",n)

mm=100
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=200
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

mm=400
v0=rand(Float64,n^2)
@btime begin ilan(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end
@btime begin tiar(Float64,nep;Neig=200,displaylevel=0,maxit=mm,tol=eps()*100,check_error_every=Inf,v=v0) end

1
