using NonlinearEigenproblems, Random, SparseArrays, PyPlot, LinearAlgebra

println("Loading the:", `gun`,"problem")
nep = nep_gallery("nlevp_native_gun")
println("Shifting and scaling the problem")
include("shift_and_scale_gun.jl");


# run tiar
println("Run TIAR")
#λ1,_,err_iar=tiar(Dnep;neigs=Inf,logger=1,maxit=30,tol=1e-5,check_error_every=1,errmeasure=err_measure)
# run Infinite Lanczos
#println("Run the infinite Lanczos method")


# RUN ILAN (IAR)
λ,W,err,_=ilan(Dnep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=1,errmeasure=err_measure,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))


pygui(true)
m,p=size(err);
for j=1:p sort!(view(err,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
ylim(ymax=10)

# plot the spectrum
#θ=range(0,stop=2*π,length=100); r=50000; c=250^2; plot(c.+r*cos.(θ),r*sin.(θ),label="region of interest")
#plot(real(λ0.+α*λ1),imag(λ0.+α*λ1),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="TIAR")
#plot(real(λ0.+α*λ;),imag(λ0.+α*λ;),marker="o",markerfacecolor=:none,c=:green,linestyle=:none,label="INF. LAN.")
#legend()

#println("Number of computed eigenpairs: ", length(λ))
#for j=1:length(λ)
#    println("Residual of the eigepair ", j, "th = ",err_measure(λ[j],W[:,j]))
#end
