using NonlinearEigenproblems,LinearAlgebra
#a=1;
#c=299792458;
#cel      = a_lat/(2*pi);
#eta=2*pi*c/a;

#gam_1=0.005*eta; # according to text?x


eps_oo_1=1.0;
om_d_1=1.1;
gam_1=0.05; # According to config file


# Normalized units so that 2*pi*c/a=1
a_lat=50;

cel      = a_lat/(2*pi);
epsf     = 8.854187817e-3;
muf      = 400*pi;
nm       = 2*pi/(a_lat*sqrt(epsf*muf));
epsilon0 = epsf*nm;
mu0      = muf*nm;
nrm     = a_lat/(2*pi*cel);

# Nrmalize
#d_sq           = d_sq          * a_lat;
#space2pml      = space2pml     * a_lat;
#pmlsize        = pmlsize       * a_lat;
#kx             = kx            * Pi/a_lat;
#eig_target_re  = eig_target_re / nrm;
#eig_target_im  = eig_target_im / nrm;
#eig_min_re     = eig_min_re    / nrm;
#eig_max_re     = eig_max_re    / nrm;
#eig_min_im     = eig_min_im    / nrm;
#eig_max_im     = eig_max_im    / nrm;
om_d_1         = om_d_1        / nrm;
gam_1          = gam_1         / nrm;



f1=s-> one(s)
f2=s-> (s-gam_1*one(s))\(-eps_oo_1*s^3+gam_1*eps_oo_1*s^2-om_d_1^2*s)
f3=s-> -s^2


include("/home/jarl/jobb/src/neppack_master2/NonlinearEigenproblems.jl/src/gallery_extra/petsc_naive_bin_read.jl")
pt="/home/jarl/archive/onelab/models-master/"

using MAT;
#A3=matread(pt*"/NonLinearEVP/mat_M15.mat")["mat_M15"]
#A2=matread(pt*"/NonLinearEVP/mat_M16.mat")["mat_M16"]
#A1=matread(pt*"/NonLinearEVP/mat_M17.mat")["mat_M17"]
A3=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M15.m.bin")
A2=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M16.m.bin")
A1=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M17.m.bin")

#A1=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M11.m.bin")
#A2=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M12.m.bin")
#A3=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M13.m.bin")

nep=SPMF_NEP([A1,A2,A3], [f1,f2,f3]);

#s0=0.259-0.00777485im;
#s0=0.259im+0.00777
s0=0.00777+0.25924im

#s0=0.00777485+0.259235im
n=size(nep,1);
v0=ones(n);
#include("/tmp/eigvec.jl");
norm(compute_Mlincomb(nep,s0,x))
#v0=x;
v0=normalize(compute_Mder(nep,s0)\v0);
#v0=normalize(compute_Mder(nep,s0)\v0);
#v0=normalize(compute_Mder(nep,s0)\v0);
#v0=normalize(compute_Mder(nep,s0)\v0);
##(s,v)=quasinewton(nep,λ=s0,v=v0,logger=1,maxit=1000)
#
(s,v)=quasinewton(nep,λ=s0,v=v0,logger=1,maxit=1000,tol=1e-9,armijo_factor=0.5)



#O=float.([0.01 ; 1; 1+1im ;1im]);
#nleigs(Float64,nep,O,logger=1)

#
#B0=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M11.m.bin")
#B1=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M12.m.bin")
#B2=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M13.m.bin")
#B3=naive_petsc_read(pt*"/NonLinearEVP/file_mat_M14.m.bin")
#
#
##s0=0.00193023+0.374668im
#s0=0.259im+0.0077
#
#pep=PEP([B0,B1,B2,B3]);
#n=size(pep,1);
#v0=ones(n);
#(s,v)=quasinewton(pep,λ=s0,v=v0,logger=1,maxit=1000,tol=1e-12)
#
