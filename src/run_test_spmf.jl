#  This is the first code in NEP-pack
# Verify the code for SPMF= sum of products of matrices and functions
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes

n=5;
A0=randn(5,5);
A1=randn(5,5);

nep0=DEP([A0,A1])

minusop= S-> -S
oneop= S -> eye(S)
expmop= S -> expm(full(-S))
fi=[minusop, oneop, expmop];

nep1=SPMF_NEP([eye(n),A0,A1],fi)


n1=norm(compute_Mder(nep0,10+1im)-compute_Mder(nep1,10+1im))
println("Test 1:",n1)

S=randn(3,3);
V=randn(n,3);
n2=norm(compute_MM(nep0,S,V)-compute_MM(nep1,S,V))
println("Test 2:",n2)

