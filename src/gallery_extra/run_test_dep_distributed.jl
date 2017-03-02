#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

nep=nep_gallery("dep_distributed");


N=5
xv,wv=gauss_legendre_weights(N,-1,0);
F=sum(wv.*sin(pi*xv)); 
println("Integration of sin(pi*x) Gauss-Legendre with N=",N,", Error:",  F- (-2/pi))

N=10
xv,wv=gauss_legendre_weights(N,-1,0);
F=sum(wv.*sin(pi*xv)); 
println("Integration of sin(pi*x) Gauss-Legendre with N=",N,", Error:",  F- (-2/pi))

N=30
xv,wv=gauss_legendre_weights(N,-1,0);
F=sum(wv.*sin(pi*xv)); 
println("Integration of sin(pi*x) Gauss-Legendre with N=",N,", Error:",  F- (-2/pi))


@time λ,v= aug_newton(nep,λ=complex(3.0),v=ones(3),displaylevel=1);

@time λ_iar,v_iar=iar(nep,displaylevel=1,maxit=20,σ=complex(3.0));

exact=2.726146249832675

println("reference error aug_newton:", abs(λ-exact))
println("reference error iar:", minimum(abs(λ_iar-exact)))


