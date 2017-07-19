workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using gplot_module

# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

nep=nep_gallery("dep0",100)
#nep=nep_gallery("pep0");


compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)


m=100; p=3; Neig=50;
try
    λ,Q,err = iar(nep,maxit=m,Neig=Neig,σ=2.0,γ=3,displaylevel=1,p=p);
    errormeasure=default_errmeasure(nep);
    for i=1:length(λ)
        println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
    end

    gcloseall()
    for i=1:m
        gsemilogy(p:p:m, err[p:p:m,i])
    end
    show_gplot()
catch e
    println(typeof(e))
    λ=e.λ
    err=e.errmeasure
    m=length(λ)
    println("Current approximations")
    for j=1:m
        println("Eigenvalue ", λ[j]," with error ", err[j])
    end
end

λ,Q,err = iar(nep,maxit=m,Neig=3,σ=0,γ=3,displaylevel=1,p=p);
try
    λ,Q,err = iar(nep,maxit=m,Neig=Neig,σ=λ[1],γ=3,displaylevel=1,p=p);
catch e
    println(e.msg)
end
