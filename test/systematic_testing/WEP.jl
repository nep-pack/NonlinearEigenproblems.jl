using GalleryWaveguide;
name="WEP";
nep=nep_gallery(WEP,nx=3,nz=3);
(λ,v)=augnewton(nep,λ=-2.8-3.1im,v=ones(size(nep,1)),tol=1e-14)
v=v/norm(v);
vapprox=v+0.03*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.1*λ;

algtol=eps()*500;
testtol=1e-10;
