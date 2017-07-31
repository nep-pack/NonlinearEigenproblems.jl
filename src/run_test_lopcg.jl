workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot
using PyCall
import Base.*



# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

nep=nep_gallery("dep0",10)
#nep=nep_gallery("pep0");


compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)


λ,Q,err = iar(nep,maxit=100,Neig=3,σ=0,γ=1,displaylevel=1,check_error_every=3);

errormeasure=default_errmeasure(nep);

AA=compute_Mder(nep,λ[1]);
nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv)
A=v->compute_Mlincomb(nept,λ[1],compute_Mlincomb(nep,λ[1],v))


n=size(nep,1)
x=u=rand(n)
q=sqrt(x⋅u)
x=x/q;  u=u/q;
v=A(x);
ρ=x⋅v;
p=zeros(n,1)

k=1
max_it=100
res=zeros(max_it)
while k<max_it
  g=v-ρ*u

  aa=[x -g p]'*[v -A(g) A(p)]; aa=(aa+aa')/2;
  mm=[x -g p]'*[u -g p]; mm=(mm+mm')/2;

  D,V=eig(aa,mm)
  ii=indmin(abs(D));
  ρ=D[ii]; δ=V[:,ii];

  p=[-g p]*δ[2:end];

  x=u=δ[1]*x + p;
  q=sqrt(x'*u);
  x/=q; u/=q;
  v=A(x);
  k+=1

  res[k]=errormeasure(λ[1],x);
end


loglog(1:max_it, res, color="red", linewidth=2.0, linestyle="--")

function *(F::Function,v::AbstractVector) F(v)
