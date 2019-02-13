# Setup a galleryproblem
name="qdep1";
nep=nep_gallery(name);
rng = MersenneTwister(1234);
(λ,v)=augnewton(nep,λ=1+1im,v=ones(size(nep,1)),tol=1e-10);
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.001*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ-0.05*λ;

algtol=eps()*100;
testtol=eps()*100;
