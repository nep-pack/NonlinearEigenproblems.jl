push!(LOAD_PATH, pwd())	# looks for modules in the current directory

using PyPlot
using PyCall

using LinearAlgebra

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver


n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];

mm=80;  # number of iterations

function myexpm(A::Array{T,2}) where {T<:Number}
    println("call expm with ",typeof(A),"\n");

    A=Array{ComplexF64,2}(A);
    F=zeros(T,size(A,1),size(A,2))
    if (size(A)==(1,1))
        F[:]=exp(A[1,1]);
        return F
    end
    Bi=eye(T,size(A,1),size(A,2))
    for k=0:100
        F=F+Bi/factorial(real(T(k)));
        Bi=Bi*A;
    end
    #F=Array{ComplexF64,2}(F);
    err=norm(expm(A)-F,1)/norm(F,1);
    if(err>eps()*100)
        println("Warning: error large:",err, " size:",size(A), " norm(A):",norm(A));

    end

#    F=expm(A);
    return F
end


nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->eye(λ),λ->myexpm(-λ)])
#nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])

function compute_y0(x,y,nep,a,b)
   T=(n,x)->cos(n*acos(x));
   U=(n,x)->n+1;
   n,N=size(x);
   y0=zeros(n,1);
   A0=nep.A[2];
   A1=nep.A[3];
   # hardcoded for 2dep
   τ=1;

   y0=A0*sum(y,2);
   for i=1:N-1
      y0=y0+(2*i/a)*U(i-1,1)*x[:,i+1];
   end

   for i=1:N+1
      y0=y0+T(i-1,1+2*τ/a)*A1*y[:,i];
   end

   y0=-(A0+A1)\y0;
end

v0=randn(n);
λ,Q,err,V,H = iar_chebyshev(nep, maxit=mm,Neig=20,σ=0.0,γ=1,displaylevel=1,check_error_every=1,compute_y0=compute_y0,v=v0);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
clf();
for i=1:m
   semilogy(1:1:m, err[1:1:m,i],color="black")
end
ylim(1e-16,1e1)




λ,Q,err,V,H = iar_chebyshev(nep, maxit=mm,Neig=20,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
   semilogy(1:1:m, err[1:1:m,i],color="red")
end
ylim(1e-16,1e1)
