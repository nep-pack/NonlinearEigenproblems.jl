workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
include("method_cork.jl");
using NEPSolver
using NEPSolver_CORK
using NEPCore



target=TargetType(0,0);
n=10;

n=20

A0=ones(n,n)
A1=eye(n,n)
A2=eye(n,n);
A3=eye(n,n)
B0=eye(n,n)
B0=B0[:,end:-1:1]
B1=tril(ones(n,n))
B2=triu(ones(n,n))

d=4

si=ones(d-1);
be=ones(d-1);
xi=ones(d-1);
M = cat(2,diagm(0 => si), zeros(d-1,1)) + cat(2,zeros(d-1,1), diagm(0 => be))
N = eye(d-1,d) + cat(2,zeros(d-1), diagm(0 => be./xi))


#M=eye(d-1,d);
#N=eye(d-1,d)*3;
L=CORKInputType([A0,A1,A2,A3],[B0,B1,B2,B2],M,N)
#
#
#
cork(L,4,target)
#
3
