using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise, CSV, DelimitedFiles
using Pkg
Pkg.activate(".")

println("Loading the:", `gun`,"problem")
nep = nep_gallery("nlevp_native_gun")
println("Shifting and scaling the problem")
include("shift_and_scale_gun.jl");


# COMPUTE REFERENCE EIGENVALUES
#println("COMPUTE REFERENCE EIGENVALUES")
λ,_,err_iar=tiar(Dnep;neigs=Inf,logger=1,maxit=200,tol=1e-5,check_error_every=Inf,errmeasure=err_measure)
CSV.write("ILAN_figures/figure_2/gun_eigs.csv", λ0.+α*λ)


# RUN ILAN (IAR)
v0=ones(size(nep,1));
λ1,W,err1,_=ilan(Dnep,σ=0,γ=1;v=v0,neigs=100,logger=1,maxit=100,tol=1e-8,check_error_every=1,errmeasure=err_measure,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
CSV.write("ILAN_figures/figure_2/gun_ilan_iar50.csv", λ0.+α*λ1)
m,p=size(err1);
for i=1:size(err1,1) for j=1:size(err1,2)	if err1[i,j]==1 err1[i,j]=NaN end end end
for j=1:p sort!(view(err1,1:m,j);rev=true) end
writedlm( "ILAN_figures/figure_2/gun_ilan_iar50_err.csv", [1:m err1], ',')


# RUN ILAN (IAR)
λ1,W,err1,_=ilan(Dnep,σ=0,γ=1;v=v0,neigs=100,logger=1,maxit=100,tol=1e-8,check_error_every=1,errmeasure=err_measure,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=200))
CSV.write("ILAN_figures/figure_2/gun_ilan_iar200.csv", λ0.+α*λ1)
m,p=size(err1);
for i=1:size(err1,1) for j=1:size(err1,2)	if err1[i,j]==1 err1[i,j]=NaN end end end
for j=1:p sort!(view(err1,1:m,j);rev=true) end
writedlm( "ILAN_figures/figure_2/gun_ilan_iar200_err.csv", [1:m err1], ',')
# increase maxit to 100 and save

# pygui(true)
# m,p=size(err);
# for j=1:p sort!(view(err,1:m,j);rev=true) end
# for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
# ylim(ymax=10)

# plot the spectrum
#θ=range(0,stop=2*π,length=100); r=50000; c=250^2; plot(c.+r*cos.(θ),r*sin.(θ),label="region of interest")
#plot(real(λ0.+α*λ1),imag(λ0.+α*λ1),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="TIAR")
#plot(real(λ0.+α*λ;),imag(λ0.+α*λ;),marker="o",markerfacecolor=:none,c=:green,linestyle=:none,label="INF. LAN.")
#legend()

#println("Number of computed eigenpairs: ", length(λ))
#for j=1:length(λ)
#    println("Residual of the eigepair ", j, "th = ",err_measure(λ[j],W[:,j]))
#end
