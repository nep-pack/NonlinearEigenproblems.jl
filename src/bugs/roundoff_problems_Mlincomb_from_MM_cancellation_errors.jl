# Illustrating scaling difficulties with compute_Mlincomb_from_MM
workspace();
push!(LOAD_PATH, pwd())
push!(LOAD_PATH, ".." )
using NEPCore
using NEPTypes

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

srand(0);
# Setup the problem
n=100
A=randn(n,n)
B=randn(n,n)
tau=2.0
dep=DEP([A,B],[0,tau])


## Compute lin comb for small m-value


# Setup a the coeff vector
m=5;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

betav=ones(alphav);
x1=compute_Mlincomb_from_MM!(dep,0,V,alphav+betav);
x2=compute_Mlincomb_from_MM!(dep,0,V,betav);
x=x1-x2;

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))


# Setup a the coeff vector
m=10;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

betav=ones(alphav);
x1=compute_Mlincomb_from_MM!(dep,0,V,alphav+betav);
x2=compute_Mlincomb_from_MM!(dep,0,V,betav);
x=x1-x2;

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))





# Setup a the coeff vector
m=20;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

betav=ones(alphav);
x1=compute_Mlincomb_from_MM!(dep,0,V,alphav+betav);
x2=compute_Mlincomb_from_MM!(dep,0,V,betav);
x=x1-x2;

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))




# Setup a the coeff vector
m=50;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

betav=ones(alphav);
x1=compute_Mlincomb_from_MM!(dep,0,V,alphav+betav);
x2=compute_Mlincomb_from_MM!(dep,0,V,betav);
x=x1-x2;

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))







# Setup a the coeff vector
m=100;
alphav=(0.6).^(1:m)
X=randn(n,m)
V,R=qr(X)

# Compute with default method
xx=DEP_Mlincomb_high_precision(A,B,tau,V,alphav)

betav=ones(alphav);
x1=compute_Mlincomb_from_MM!(dep,0,V,alphav+betav);
x2=compute_Mlincomb_from_MM!(dep,0,V,betav);
x=x1-x2;

norm(x-xx)/norm(xx)
println("m=",m)
println("Error:",Float64(norm(x-xx)/norm(xx)))
