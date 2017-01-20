#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery


println("Load dep0")
nep=nep_gallery("dep0")


################# NEWTON #########################3


println("Running Newton on random dep")
nep=nep_gallery("dep0")

λ=NaN;
x=NaN
try
    λ,x =newton(nep,displaylevel=1);
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(λ)
println(compute_resnorm(nep,λ,x))




################# INFINITE ARNOLDI #####################
n=nep.n;
m=10;
V=zeros(n*(m+1),m+1);
H=zeros(m+1,m);
Bv=zeros(n,m+1);

V[1:n,1]=rand(n,1);
V[1:n,1]=V[1:n,1]/norm(V[1:n,1]);

alpha=zeros(m+1);

Minv=LinSolver(compute_Mder(nep,0.0));

for k=1:m
 # action of BB
 W=zeros(n,k+1);
 W[:,2:k+1]=reshape(V[1:n*k,k],n,k);
 
	
 alpha[1]=0; alpha[2:k+1]=1;
 Bv[1:n,1]=compute_Mlincomb_from_Mder(nep,0.0,W,alpha[1:k+1]);
 Bv[1:n,1]=-Minv.solve(Bv[1:n,1]);
 for j=2:k+1
  Bv[:,j]=W[:,j]/j;
 end
 
 # new vector
 vv=reshape(Bv[:,1:k+1],(k+1)*n,1);

 # double GS-orth 
 h=V[1:(k+1)*n,1:k]'*vv;

 vv=vv-V[1:(k+1)*n,1:k]*h;
 
 g=V[1:(k+1)*n,1:k]'*vv;
 vv=vv-V[1:(k+1)*n,1:k]*g;

 H[1:k,k]=h+g;


 beta=norm(vv);

 H[k+1,k]=beta;


 V[1:(k+1)*n,k+1]=vv/beta;


end


D,V=eig(H[1:m,1:m]);
#D=D[1];
D=1./D;



err=abs(λ-D);
minimum(err);


