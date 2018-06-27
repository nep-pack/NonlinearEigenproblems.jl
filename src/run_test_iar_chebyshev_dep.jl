workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory

using PyPlot
using PyCall

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery



#explicit import needed for overloading functions from packages
import NEPCore.compute_Mlincomb


nep=nep_gallery("dep0_tridiag",1000)

compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)


function compute_y0(x,y,nep,a,b)
   T=(n,x)->cos(n*acos(x));
   c=(a+b)/(a-b);   k=2/(b-a);

   N=size(x,2);
   n=size(x,1);


   y0=zeros(n,1);
   for i=1:N
       y0=y0+T(i-1,c)*x[:,i];
   end

   for j=1:length(nep.A)
       for i=1:N+1
           y0=y0-T(i-1,-k*nep.tauv[j]+c)*(nep.A[j]*y[:,i]);
       end
   end

   y0=sum(nep.A)\y0;

end

(λ,Q,err,V,H)=iar_chebyshev(nep,σ=0,γ=1,Neig=30,v=ones(4),displaylevel=0,maxit=100,tol=eps()*100,check_error_every=1,compute_y0=compute_y0)

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(3:3:m, err[3:3:m,i],color="black")
end
ylim(1e-16,1e1)
