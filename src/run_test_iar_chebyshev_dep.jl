workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot
using PyCall

n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];

#a=-1; b=0; k=2/(b-a); c=(a+b)/(a-b);
#cheb_vect=(t,Z)->cos((0:(size(Z,2)-1))*acos(t))';
#cheb2_vect_m1=(Z)->(0:(size(Z,2)-1))';
#Mterm=(t,X)->k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));
#y0comp=(X,Y)->(A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));

# Chebyshev polynomials of the first kind
function T(i,x)
   if i==0
      1
   elseif i==1
      x
   else
      2*x*T(i-1,x)-T(i-2,x)
   end
end

# Chebyshev polynomials of the second kind
function U(i,x)
   if i==0
      1
   elseif i==1
      2*x
   else
      2*x*U(i-1,x)-U(i-2,x)
   end
end

# define the quadratic DEP
# IDEA: should we define PDEP=polynomial DEP?
nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->1,λ->exp(-λ)])

function compute_y0(x,y,nep)
   n,N=size(y);
   y0=zeros(n,1);
   A0=nep.A[2];
   A1=nep.A[3];


   for i=1:N
      y0=y0+T(i,3)*y[:,i];
   end
   y0=-A1*y0;

   for i=1:N-1
      # be careful with the index of x. It starts from zero.
      y0=y0+2*i*U(i-1,1)*x[:,i+1];
   end
   y0=y0-A0*sum(y,2);
   y0=(A0+A1)\y0;

end



λ,Q,err = iar_chebyshev(nep, compute_y0, maxit=100,Neig=10,σ=2.0,γ=3,displaylevel=1,check_error_every=3);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(3:3:m, err[3:3:m,i],color="black")
end
