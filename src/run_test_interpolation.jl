#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery


println("Running Newton on random dep")
nep=nep_gallery("dep0")

λ1=NaN;
x1=NaN
try
    λ1,x1 =newton(nep,displaylevel=1,maxit=40, λ=-0.75, v=ones(size(nep,1),1));
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ1=e.λ
    x1=e.v
end
println("Resnorm: ",compute_resnorm(nep,λ1,x1))
println(λ1)
println(x1)


λ2=NaN;
x2=NaN
println("\nRunning Newton on interpolated dep")
intpoints = [λ1-1, λ1, λ1+1.5]
pep = interpolate(nep, intpoints)
try
    λ2,x2 =newton(pep,displaylevel=1,maxit=40, λ=-0.75, v=ones(size(nep,1),1));
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ2=e.λ
    x2=e.v
end
println("Resnorm: ",compute_resnorm(pep,λ2,x2))
println(λ2)
println(x2)

println("\nDifferences\nEigenvalue ", abs(λ1-λ2))
println("Error norm = ",norm(x1-x2), " Eigenvector norm = ", norm(x1))

println("\nRunning Newton on interpolated dep (higher degree of interpolating pep)")
#intpoints = [λ1-7, λ1-5, λ1-2, λ1-1, λ1, λ1+1, λ1+2, λ1+5, λ1+7]
intpoints = [λ1-5, λ1-1, λ1, λ1+5, λ1+1, λ1+5im, λ1+1im, λ1-5im, λ1-1im]
pep = interpolate(nep, intpoints)
try
    λ2,x2 =newton(pep,displaylevel=1,maxit=40, λ=-0.75, v=ones(size(nep,1),1));
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ2=e.λ
    x2=e.v
end
println("Resnorm: ",compute_resnorm(pep,λ2,x2))
println(λ2)
println(x2)

println("\nDifferences\nEigenvalue ", abs(λ1-λ2))
println("Error norm = ",norm(x1-x2), " Eigenvector norm = ", norm(x1))


println("\n\n\n\n-------------------------------------------")
println("\n\n\n\nRunning Newton on random pep (original pep of degree 2)")
nep=nep_gallery("pep0")
try
    λ1,x1 =newton(nep,displaylevel=1,maxit=40, λ=1, v=ones(size(nep,1),1));
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ1=e.λ
    x1=e.v
end
println("Resnorm: ",compute_resnorm(nep,λ1,x1))
println(λ1)
#println(x1)

println("\nRunning Newton on interpolated pep (interpolation of degree 2)")
intpoints = [λ1-1, λ1, λ1+1]

pep = interpolate(nep, intpoints)
try
    λ2,x2 =newton(pep,displaylevel=1,maxit=40, λ=1, v=ones(size(nep,1),1));
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ2=e.λ
    x2=e.v
end
println("Resnorm: ",compute_resnorm(pep,λ2,x2))
println(λ2)
#println(x2)
println("\nCoefficient matrix differences (monomes): ")
println(norm(nep.A[1]-pep.A[1]))
println(norm(nep.A[2]-pep.A[2]))
println(norm(nep.A[3]-pep.A[3]))
#println(nep.A[1])
#println(pep.A[1])


println("\nDifferences\nEigenvalue ", abs(λ1-λ2))
println("Error norm = ",norm(x1-x2), " Eigenvector norm = ", norm(x1))



println("\nRunning Newton on interpolated pep (interpolation of degree 4)")
intpoints = [λ1-3, λ1-1, λ1, λ1+1, λ1+3]
pep = interpolate(nep, intpoints)
try
    λ2,x2 =newton(pep,displaylevel=1,maxit=40, λ=1, v=ones(size(nep,1),1));
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ2=e.λ
    x2=e.v
end
println("Resnorm: ",compute_resnorm(pep,λ2,x2))
println(λ2)
#println(x2)
println("\nCoefficient matrix differences (monomes): ")
println(norm(nep.A[1]-pep.A[1]))
println(norm(nep.A[2]-pep.A[2]))
println(norm(nep.A[3]-pep.A[3]))
println(norm(pep.A[4]))
println(norm(pep.A[5]))


println("\nDifferences\nEigenvalue ", abs(λ1-λ2))
println("Error norm = ",norm(x1-x2), " Eigenvector norm = ", norm(x1))




