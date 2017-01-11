#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore

println("Test Newton")

n=5;
srand(0) # reset the random seed
A0=randn(n,n);
A1=randn(n,n);
I=eye(n,n);
tau=1;

#nep.n=n;
function DEP_Md(λ,i=0)
    # short-hand definition of functions
    if (i==0)
       return -λ*I+A0+A1*exp(-tau*λ)
    else
       return -I-tau*A1*exp(-tau*λ)
    end
end

nep=NEP(n,DEP_Md);


λ=NaN;
x=NaN

    v = aug_newton(nep);
    v









