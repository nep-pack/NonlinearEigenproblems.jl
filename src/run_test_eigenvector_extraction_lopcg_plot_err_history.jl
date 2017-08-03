workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

nep=nep_gallery("dep0_sparse",100)
nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv)

#nep=nep_gallery("pep0",50);
#nept=PEP([nep.A[1]',nep.A[2]',nep.A[3]'])

λ,Q,err = iar(nep,maxit=100,Neig=2,σ=1.0,γ=1,displaylevel=1,check_error_every=1);

errormeasure=default_errmeasure(nep);

A=v->compute_Mlincomb(nept,conj(λ[1]),compute_Mlincomb(nep,λ[1],v))

n=size(nep,1)

#x=rand(Complex128,n)
x=Q[:,1]+1e-6*rand(Complex128,n);
q=zeros(Complex128,n)

x/=norm(x);
v=A(x);
ρ=x⋅v;


k=1
max_it=convert(Int,1e5)
res_vec=zeros(Float64,max_it)
err=1
tol=1e-12
while (k<max_it)&&(err>tol)
  println("Iteration: ",k/max_it)
  g=v-ρ*x

  aa=[x -g q]'*[v -A(g) A(q)]; aa=(aa+aa')/2;
  mm=[x -g q]'*[x -g q]; mm=(mm+mm')/2;

  D,V=eig(aa,mm)
  ii=indmin(abs(D));
  ρ=D[ii]; δ=V[:,ii];

  q=[-g q]*δ[2:end];

  x=δ[1]*x+q;
  x/=norm(x);
  v=A(x);
  k+=1

  err=errormeasure(λ[1],x)
  res_vec[k]=errormeasure(λ[1],x);

end

using PyPlot
using PyCall
loglog(1:max_it, res_vec, color="red", linewidth=2.0, linestyle="--")
