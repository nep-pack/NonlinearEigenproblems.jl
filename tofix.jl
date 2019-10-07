using NonlinearEigenproblems, LinearAlgebra
nep=nep_gallery("dep_symm_double",10);
v0=ones(size(nep,1));
λ,v=ilan(nep;v=v0,tol=1e-5,neigs=3);
norm(compute_Mlincomb!(nep,λ[1],v[:,1])) # Is it an eigenvalue?
λ    # print the computed eigenvalues
