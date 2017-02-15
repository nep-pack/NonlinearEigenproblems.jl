#  This is the first code in NEP-pack
# Verify the code for SPMF= sum of products of matrices and functions
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes

n=5;
srand(1)
A0=randn(5,5);
A1=randn(5,5);
t=3.0
nep0=DEP([A0,A1],[0,t])

minusop= S-> -S
oneop= S -> eye(S)
expmop= S -> expm(full(-t*S))
fi=[minusop, oneop, expmop];

nep1=SPMF_NEP([eye(n),A0,A1],fi)


n1=norm(compute_Mder(nep0,10+1im)-compute_Mder(nep1,10+1im))
println("Test 1:",n1)

S=randn(3,3);
V=randn(n,3);
n2=norm(compute_MM(nep0,S,V)-compute_MM(nep1,S,V))
println("Test 2:",n2)


v0=randn(n)
println("Resinv test 0")
res_inv(nep0,displaylevel=1,v=v0)
println("Resinv test 1")
res_inv(nep1,displaylevel=1,v=v0)


n3=norm(compute_Mder(nep0,-1,5)-compute_Mder_from_MM(nep1,-1,5))
println("Test 3: comparing nep0 with nep1 norm=",n3)
λ=1+0.3im
n4=norm(compute_Mder(nep0,λ,20)-compute_Mder_from_MM(nep0,λ,20))
println("Test 4: comparing Mder with Mder_from_MM nep0. norm=",n4)


n5=norm(compute_Mder(nep0,λ,0)-compute_Mder(nep1,λ,0))
println("Test 5:  comparing Mder (i=0) nep0 and nep1. norm=",n5)

n6=norm(compute_Mder(nep0,λ,1)-compute_Mder(nep1,λ,1))
println("Test 6:  comparing Mder (i=1) nep0 and nep1. norm=",n6)

n7=norm(compute_Mder(nep0,λ,15)-compute_Mder(nep1,λ,15))
println("Test 7:  comparing Mder (i=15) nep0 and nep1. norm=",n7)
