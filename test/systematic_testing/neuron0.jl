# Setup a galleryproblem
name="neuron0";
nep=nep_gallery(name);
rng = MersenneTwister(1234);
(λ,v)=augnewton(Float64,nep,λ=-4,v=ones(size(nep,1)));
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.001*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ-0.05*λ;

algtol=1e-10;
testtol=1e-10;
