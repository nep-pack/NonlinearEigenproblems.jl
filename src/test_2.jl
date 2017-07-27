workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
#using gplot_module


nep=nep_gallery("pep0",3)
#nep=nep_gallery("dep0",3)

#a,b,c,d=rand(4,1)+rand(4,1)*im
#a,b,c,d=[1,0,0,1]
γ=-1+im;
a,b,c,d=[1,-γ,1,conj(γ)]    # Cayley
#a,b,c,d=[10,-1,1,10]



nep1=MobiusTransformNEP(nep,a=a,b=b,c=c,d=d);

compute_Mlincomb(nep::MobiusTransformNEP,λ::Number,V;a=ones(size(V,2)))= compute_Mlincomb_from_MM!(nep,λ,V,a)
k=100; n=size(nep,1);
σ=0; α=1.^(0:k); α[1]=0; y=rand(n,k+1);
z=compute_Mlincomb(nep1,σ,y[:,1:k+1],a=α[1:k+1]);
z
