using NonlinearEigenproblems, Random, SparseArrays, PyPlot, LinearAlgebra

println("Loading the:", `gun`,"problem")
nep = nep_gallery("nlevp_native_gun")
println("Shifting and scaling the problem")
include("shift_and_scale_gun.jl");


# run tiar
println("Run TIAR")
λ1,_,err_iar=tiar(Dnep;Neig=200,displaylevel=0,maxit=100,tol=1e-5,check_error_every=1,errmeasure=err_measure)
# run Infinite Lanczos
println("Run the infinite Lanczos method")
λ,W=ilan(Dnep,Neig=80,displaylevel=0,maxit=mm,tol=1e-6,check_error_every=Inf,errmeasure=err_measure)

# plot the spectrum
θ=range(0,stop=2*π,length=100); r=50000; c=250^2; plot(c.+r*cos.(θ),r*sin.(θ),label="region of interest")
plot(real(λ0.+α*λ1),imag(λ0.+α*λ1),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="TIAR")
plot(real(λ0.+α*λ;),imag(λ0.+α*λ;),marker="o",markerfacecolor=:none,c=:green,linestyle=:none,label="INF. LAN.")
legend()

println("Number of computed eigenpairs: ", length(λ))
for j=1:length(λ)
    println("Residual of the eigepair ", j, "th = ",err_measure(λ[j],W[:,j]))
end
