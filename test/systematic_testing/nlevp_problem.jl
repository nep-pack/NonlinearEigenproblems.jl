using GalleryNLEVP;
nep=nep_gallery(NLEVP_NEP,name)
if (tol0<0)
   tol0=1e-10;
end
@show name
(λ,v)=augnewton(nep,λ=λ0,v=ones(size(nep,1)),displaylevel=1,tol=tol0);
@show λ
vapprox=v+0.001*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ-0.05*λ;
if algtol<0;
    algtol=1e-10;
end
if algtol<0
  testtol=1e-10;
end
algtol=-1;
testtol=-1;
tol0=-1;
