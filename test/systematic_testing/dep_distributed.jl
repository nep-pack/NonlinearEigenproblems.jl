name="dep_distributed";
nep=nep_gallery(name);
(λ,v)=augnewton(nep,λ=-0.4+1im,v=ones(size(nep,1)))
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.01*λ;
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                eps()*100,1e-10,"dep_distributed");
push!(tests,neptest)
