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

nep=nep_gallery("dep0",1000)
#nep=nep_gallery("pep0");


compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)


try
    λ,Q,err = iar(nep,maxit=100,Neig=50,σ=2.0,γ=3,displaylevel=1,check_error_every=3);
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
    println(e.msg)
    println("Current approximations")
    for j=1:m
        println("Eigenvalue ", λ[j]," with error ", err[j])
    end
    global σ=λ[1]   # make this variable visible outisde catch
end

try
    λ,Q,err = iar(nep,maxit=100,Neig=3,σ=σ,γ=3,displaylevel=1,check_error_every=3);
catch e
    println(e.msg)
end
