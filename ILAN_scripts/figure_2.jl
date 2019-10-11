using NonlinearEigenproblems, Random, SparseArrays, LinearAlgebra, PyPlot

pygui(true)

println("Loading the:", `gun`,"problem")
nep = nep_gallery("nlevp_native_gun")
println("Shifting and scaling the problem")
include("shift_and_scale_gun.jl");

# COMPUTE REFERENCE EIGENVALUES
println("COMPUTE REFERENCE EIGENVALUES")
λ,_,err_iar=tiar(Dnep;neigs=Inf,logger=1,maxit=200,tol=1e-5,check_error_every=Inf,errmeasure=err_measure)


# RUN ILAN (IAR 50)
v0=ones(size(nep,1));
λ1,W,err1,_=ilan(Dnep,σ=0,γ=1;v=v0,neigs=Inf,logger=1,maxit=100,tol=1e-8,check_error_every=1,errmeasure=err_measure,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
m,p=size(err1);
for i=1:size(err1,1) for j=1:size(err1,2)	if err1[i,j]==1 err1[i,j]=NaN end end end
for j=1:p sort!(view(err1,1:m,j);rev=true) end


# RUN ILAN (IAR 200)
λ2,W,err2,_=ilan(Dnep,σ=0,γ=1;v=v0,neigs=Inf,logger=1,maxit=100,tol=1e-8,check_error_every=1,errmeasure=err_measure,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=200))
m,p=size(err1);
for i=1:size(err1,1) for j=1:size(err1,2)	if err1[i,j]==1 err1[i,j]=NaN end end end
for j=1:p sort!(view(err1,1:m,j);rev=true) end


# PLOT ERROR HIST
figure(1)
m,p=size(err2);
for j=1:p semilogy(1:m,err1[1:m,j],color="red",linestyle="--") end
semilogy(1:m,err1[1:m,1],color="red",linestyle="--",label="ILAN (IAR 50)")
for j=1:p semilogy(1:m,err2[1:m,j],color="black",linestyle="-") end
semilogy(1:m,err2[1:m,1],color="black",linestyle="-",label="ILAN (IAR 200)")
legend()

# plot the spectrum
figure(2)
θ=range(0,stop=2*π,length=100); r=50000; c=250^2; plot(c.+r*cos.(θ),r*sin.(θ),label="region of interest")
plot(real(λ0.+α*λ),imag(λ0.+α*λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="eigenvalues")
plot(real(λ0.+α*λ1),imag(λ0.+α*λ1),marker="s",markerfacecolor=:none,c=:black,linestyle=:none,label="ILAN (IAR 50)")
plot(real(λ0.+α*λ2;),imag(λ0.+α*λ2;),marker="o",markerfacecolor=:none,c=:green,linestyle=:none,label="ILAN (IAR 200)")
legend(numpoints = 1)	# display only one marker
