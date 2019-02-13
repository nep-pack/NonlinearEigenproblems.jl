# Setup a galleryproblem
name="dep0_sparse";
nep=nep_gallery(name);
rng = MersenneTwister(1234);
(λ,v)=augnewton(Float64,nep,λ=-4,v=ones(size(nep,1)),tol=1e-10);
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.01*λ;

algtol=eps()*100;
testtol=eps()*1000;
