workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver

#import NEPSolver.iar_chebyshev;
#include("../src/method_iar_chebyshev.jl");

#TODO: add this to the gallery
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

nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])

# The user can create his own orthogonalization function to use in IAR
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

   return y0=-(A0+A1)\y0;
end
# Then it is needed to create a type to access to this function
import NEPSolver.ComputeY0Cheb
import NEPSolver.AbstractPrecomputeData
abstract type ComputeY0Cheb_QDEP <: ComputeY0Cheb end
type PrecomputeData_QDEP <: AbstractPrecomputeData end

# And then introduce a function dispatch for this new type in order to use
# the defined orthogonalization function
import NEPSolver.precompute_data
function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},a,b,m,γ,σ)
    return PrecomputeData_QDEP()
end

import NEPSolver.compute_y0_cheb
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb_QDEP},x,y,M0inv,precomp::PrecomputeData_QDEP)
   return compute_y0(x,y,nep,-1,1)
end




#TODO: fix this function

v0=randn(n);
λ,Q,err = iar_chebyshev(nep,compute_y0_method=ComputeY0Cheb_QDEP,maxit=mm,Neig=10,σ=0.0,γ=1,displaylevel=1,check_error_every=1,v=v0);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
