# Illustrating scaling difficulties with compute_Mlincomb_from_MM
push!(LOAD_PATH, pwd())
push!(LOAD_PATH, ".." )
using NEPCore
using NEPTypes
using LinearAlgebra
using Random

# manually clc
for jj=1:20
    println("\n")
end

function DEP_Mlincomb_high_precision(A,B,tau,V,alphav)
    A=Array{BigFloat}(A);
    B=Array{BigFloat}(B);
    alphav=Array{BigFloat}(alphav);
    tau=BigFloat(tau)
    xx=zeros(BigFloat,n)
    xx=alphav[1]*(A+B)*V[:,1]
    xx+=alphav[2]*(-V[:,2]-tau*B*V[:,2])
    for i=2:(size(V,2)-1)
        xx+=B*V[:,i+1]*((-tau)^i)*alphav[i+1]
    end
    return xx
end

Random.seed!(0);
# Setup the problem
n=100
A=randn(n,n)
B=randn(n,n)
tau=2.0
dep=DEP([A,B],[0,tau])


## Compute lin comb for small m-value

# Setup a the coeff vector
m=10;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
x=compute_Mlincomb_from_MM(dep,0,V,alphav)
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))

## Compute lin comb for larger m-value

# Setup a the coeff vector
m=30;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
x=compute_Mlincomb_from_MM!(dep,0,V,alphav)
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))

# Setup a the coeff vector
m=100;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
alphav[2]=0;
x=compute_Mlincomb_from_MM!(dep,0,V,alphav)
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))



# Setup a the coeff vector
m=500;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
x=compute_Mlincomb_from_MM!(dep,0,V,alphav)
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))
