#  Tests for the projected problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
#using Winston # For plotting

peptest=true
local nep::NEP
if (peptest)
    nep=nep_gallery("pep0");
else
    
    n=5;
    srand(1)
    A0=randn(5,5);
    A1=randn(5,5);
    t=3.0

    minusop= S-> -S
    oneop= S -> eye(S)
    expmop= S -> expm(full(-t*S))
    fi=[minusop, oneop, expmop];
    
    nep=SPMF_NEP([eye(n),A0,A1],fi)
end

#

println("Running Newton Raphson")
λ,x =newton(nep,maxit=30,
            displaylevel=1,λ=1+1im,tolerance=1e-6);
#
λ_exact=λ
ev2=zeros(0)

pnep=Proj_NEP(nep);
V=randn(size(nep,1),2)
Q,R=qr(hcat(V,x)) # Make the eigenspace a part of the projection subspace
set_projectmatrices!(pnep,Q,Q);
println("Running Newton on projected problem with very good start value")
λ1,z1=newton(pnep,λ=(λ_exact+0.00001),displaylevel=1)

x1=Q*z1; x1=x1/x1[1];

# should be small since the eigenvector is in the subspace
println("Difference of solution from projected problem:", norm(x/x[1]-x1))

