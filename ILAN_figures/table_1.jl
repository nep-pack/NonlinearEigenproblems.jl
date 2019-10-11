using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, BenchmarkTools


mm=50
println("START3")

########################## n=100 ##########################
n=100;
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)
x = range(0, stop = π, length = n); h=x[2]-x[1]; LL=LL/(h^2); LL=-kron(LL,LL); A=LL
b=broadcast((x,y)->-x*sin(y-x),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])
nep=DEP([A,B],[0,1])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));
v0=ones(n^2)+ones(n^2)*1im


λ,v=iar(nep;maxit=mm,tol=1e-8,neigs=Inf,logger=0,check_error_every=Inf,v=v0,errmeasure=rel_err)
@btime λ,v=iar(nep;maxit=mm,tol=1e-8,neigs=Inf,logger=0,check_error_every=Inf,v=v0,errmeasure=rel_err)

λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))
@btime λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))

λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
@btime λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))

λ4,v4,err4,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=100))
@btime λ4,v4,err4,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=100))

println("λ=",length(λ))
println("λ2=",length(λ2))
println("λ3=",length(λ3))
println("λ4=",length(λ4))


########################## n=300 ##########################
n=300;
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)
x = range(0, stop = π, length = n); h=x[2]-x[1]; LL=LL/(h^2); LL=-kron(LL,LL); A=LL
b=broadcast((x,y)->-x*sin(y-x),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])
nep=DEP([A,B],[0,1])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));
v0=ones(n^2)


λ,v=iar(nep;maxit=mm,tol=1e-8,neigs=Inf,logger=0,check_error_every=Inf,v=v0,errmeasure=rel_err)
@btime λ,v=iar(nep;maxit=mm,tol=1e-8,neigs=Inf,logger=0,check_error_every=Inf,v=v0,errmeasure=rel_err)

λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))
@btime λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))

λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
@btime λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))

λ4,v4,err4,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=100))
@btime λ4,v4,err4,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=100))

println("λ=",length(λ))
println("λ2=",length(λ2))
println("λ3=",length(λ3))
println("λ4=",length(λ4))

########################## n=500 ##########################
n=500;
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)
x = range(0, stop = π, length = n); h=x[2]-x[1]; LL=LL/(h^2); LL=-kron(LL,LL); A=LL
b=broadcast((x,y)->-x*sin(y-x),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])
nep=DEP([A,B],[0,1])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));
v0=ones(n^2)


λ,v=iar(nep;maxit=mm,tol=1e-8,neigs=Inf,logger=0,check_error_every=Inf,v=v0,errmeasure=rel_err)
@btime λ,v=iar(nep;maxit=mm,tol=1e-8,neigs=Inf,logger=0,check_error_every=Inf,v=v0,errmeasure=rel_err)

λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))
@btime λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))

λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
@btime λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))

λ4,v4,err4,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=100))
@btime λ4,v4,err4,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=mm,tol=1e-8,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=100))

println("λ=",length(λ))
println("λ2=",length(λ2))
println("λ3=",length(λ3))
println("λ4=",length(λ4))
