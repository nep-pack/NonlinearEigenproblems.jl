# Setup a galleryproblem
name="nlevp_native_hadeler";
nep=nep_gallery(name);
rng = MersenneTwister(1234);
(λ,v)=augnewton(nep,λ=0.2,v=ones(size(nep,1)),tol=1e-13);
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.005*λ;

algtol=eps()*500;
testtol=1e-10;
