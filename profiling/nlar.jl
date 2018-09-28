#Profiler for the Nonlinear Arnoldi method implement in src/method_nlar.jl
using ..NEPSolver 
using Profile
using Test
using LinearAlgebra
using Random

Profile.clear()

TOL = 1e-10;

Random.seed!(0)

nep = nep_gallery("nlevp_native_gun")
n = size(nep,1);

#The eigenvalues of the gun problem lie within a semi-circle centred approximately around the "shift" with an approximate radius "scale"
shift = 250^2;
scale = 330^2-220^2;

#Run Quasi-Newton for initial guess obtained from the knowledge of the eigenvalue distribution
λ_ref,v_ref = quasinewton(nep, λ = shift+scale*(-0.131403+0.00759532im), v = ones(n), displaylevel = 0, tol = TOL/50, maxit = 500)

#Shift and scale the NEP(mainly to avoid round-off errors because of the large entries in the coefficient matrices)
nep1_spmf = SPMF_NEP(get_Av(nep),get_fv(nep))
nep1 = shift_and_scale(nep1_spmf,shift=shift,scale=scale);

#Start profiling
println("Running profiler on NLAR with n = 1e+7 and delay = 1.0")
Profile.init(n = 10^7, delay = 1.0);
#Run NLAR on the shifted and scaled NEP (nev set to 1 to save time.
@profile nlar(nep1, tol=TOL, λ0 = 0,maxit=500, nev = 5, R=0.01,displaylevel=1,v0=ones(n),eigval_sorter=residual_eigval_sorter,inner_solver_method=NEPSolver.IARInnerSolver,qrfact_orth=false,num_restart_ritz_vecs=3,max_subspace=50);
Profile.print(format=:flat,sortedby=:count);
