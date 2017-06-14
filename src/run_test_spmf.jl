#  This is the first code in NEP-pack
# Verify the code for SPMF= sum of products of matrices and functions
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using BenchmarkTools

n=5;
srand(1)
A0=randn(n,n);
A1=randn(n,n);
t=3.0
nep0=DEP([A0,A1],[0,t])

minusop= S-> -S
oneop= S -> eye(S)
expmop= S -> expm(full(-t*S))
fi=[minusop, oneop, expmop];

nep1=SPMF_NEP([eye(n),A0,A1],fi)

nep2=SPMF_NEP([eye(n),A0,A1],fi, true)

n1=norm(compute_Mder(nep0,10+1im)-compute_Mder(nep1,10+1im))
println("Test 1: ",n1)


S=randn(3,3);
V=randn(n,3);
n2=norm(compute_MM(nep0,S,V)-compute_MM(nep1,S,V))
println("Test 2: ",n2)
n2_1=norm(compute_MM(nep0,S,V)-compute_MM(nep2,S,V))
println("Test 2.1: ",n2_1)


v0=randn(n)
println("Resinv test 0")
resinv(nep0,displaylevel=1,v=v0)
println("Resinv test 1")
resinv(nep1,displaylevel=1,v=v0)


n3=norm(compute_Mder(nep0,-1,5)-compute_Mder_from_MM(nep1,-1,5))
println("Test 3: comparing nep0 with nep1 norm= ",n3)
λ=1+0.3im
n4=norm(compute_Mder(nep0,λ,20)-compute_Mder_from_MM(nep0,λ,20))
println("Test 4: comparing Mder with Mder_from_MM nep0. norm= ",n4)


n5=norm(compute_Mder(nep0,λ,0)-compute_Mder(nep1,λ,0))
println("Test 5: comparing Mder (i=0) nep0 and nep1. norm= ",n5)

n6=norm(compute_Mder(nep0,λ,1)-compute_Mder(nep1,λ,1))
println("Test 6: comparing Mder (i=1) nep0 and nep1. norm= ",n6)

n7=norm(compute_Mder(nep0,λ,15)-compute_Mder(nep1,λ,15))
println("Test 7: comparing Mder (i=15) nep0 and nep1. norm= ",n7)




println("")

n = 400
m = 20
srand(12)
A0=randn(n,n)
A1=randn(n,n)
A2=randn(n,n)
A3=randn(n,n)
A4=randn(n,n)
t=3.0

minusop = X-> -X
oneop = X -> eye(X)
expmop = X -> expm(-t*X)
logmop = X -> logm(X)
sqrtmop = X -> sqrtm(X)
compop = X -> logm(X)+sqrtm(X)

fi=[minusop, oneop, expmop, logmop, sqrtmop, compop];

nep3=SPMF_NEP([eye(n),A0,A1,A2,A3,A4],fi)
nep4=SPMF_NEP([eye(n),A0,A1,A2,A3,A4],fi, true) #ACTIVATE PRE-COMPUTATION OF SCHUR DECOMPOSITION

srand(123)
SS = randn(m,m)
Q,R = qr(SS)
D = diagm(10*rand(m))
SS = Q*D*Q'

VV = randn(n,m)

function without_fact(XS)
    return compute_MM(nep3,XS,VV)
end

function with_fact(XS)
    return compute_MM(nep4,XS,VV)
end

n8=norm(without_fact(SS)-with_fact(SS))
println("Test 8: comparing MM with and without Schur pre-factorization\n        (S is Hermitian PSD)\n        Difference norm = ",n8)

println("        Running without_fact...");
b1=@benchmark without_fact(SS)
show(STDOUT, "text/plain", b1) # Pretty print it

println("\n\n        Running with_fact...");
b2=@benchmark with_fact(SS)
show(STDOUT, "text/plain", b2) # Pretty print it


srand(1234)
SS_alt = randn(m,m); SS_alt = (SS_alt + eye(size(SS_alt,1))*1.05*norm(SS_alt))/4

n9=norm(without_fact(SS_alt)-with_fact(SS_alt))
println("\n")
println("\nTest 9: comparing MM with and without Schur pre-factorization\n        (S is real and eigenvalues have prositive real part)\n        Difference norm = ",n9)

println("        Running without_fact...");
b1=@benchmark without_fact(SS_alt)
show(STDOUT, "text/plain", b1) # Pretty print it

println("\n\n        Running with_fact...");
b2=@benchmark with_fact(SS_alt)
show(STDOUT, "text/plain", b2) # Pretty print it



