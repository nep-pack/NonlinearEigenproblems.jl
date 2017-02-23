#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

nep=nep_gallery("dep_distributed");
#
#S=randn(10,10);
#Z=distributed_kernel(S)
#
#println(Z)

λ,v= aug_newton(nep,λ=complex(3.0),v=ones(3));


exact=2.726146249832675

println("reference error:", abs(λ-exact))
