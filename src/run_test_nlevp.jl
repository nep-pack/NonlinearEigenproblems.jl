#  Example illustrating loading a NEP from the NLEVP toolbox
# 

workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPCore
using NEPSolver
using NEPSolver_MSLP
using Gallery

## Init
nep_name="gun"
println("Loading \"",nep_name,"\" from Berlin-Manchester collection")
nep=nlevp_gallery_import(nep_name)
println("Loading completed. Size of problem n=",nep.n)
println("Testing some basic operations on the nep")
M=compute_Mder(nep,160^2+3im)
Mp=compute_Mder(nep,160^2+3im,1)
println("It worked.")

## Run it

# Approximate eigenvalue of gun problem from Liao, Bai, Lee, Ko
λ0=((1.495 + 0.000021)*1e2)^2

v0=ones(size(nep,1))

println("Running aug newton")
λ,v=aug_newton(nep,λ=λ0,v=v0,
               displaylevel=2,maxit=30,tolerance=1e-6)
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))
println("Running newton")
λ,v=newton(nep,λ=λ0,v=v0,
               displaylevel=2,maxit=30,tolerance=1e-6)
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))

println("Running MSLP")
λ,v=mslp(nep,λ=λ0,
         displaylevel=2,maxit=30,tolerance=1e-4,eigsolver="inv_it")
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))



println("Running one step of MSLP and then aug_newton")
λ1=λ0
v1=copy(v0)
try 
    λ1,v1=mslp(nep,λ=λ0,
               displaylevel=2,maxit=1,tolerance=-1,
               eigsolver="inv_it")
catch e
    λ1=e.λ
    v1=copy(e.v)
end
v1=v1/norm(v1)
#println("sqrt(λ1)=",sqrt(λ1))
λ,v=aug_newton(nep,λ=λ1,v=v1,c=v0,
               displaylevel=2,maxit=30,tolerance=1e-6)

println("Found eigenvalue \sqrt{λ}=",sqrt(λ))




println("Running resinv")
λ,v=res_inv(nep,λ=λ0,
               displaylevel=2,maxit=20,tolerance=1e-4)
println("Found eigenvalue \sqrt{λ}=",sqrt(λ))

