"""
    Module with the Safeguarded Iteration
"""
module NEPSolver_SG_ITER

    using NEPCore
    using MATLAB   
    export sg_iteration



    #############################################################################
    # Safeguarded iteration
    function sg_iteration(nep::NEP;
                  errmeasure::Function =
                  default_errmeasure(nep::NEP),
                  tol_outer = 1e-2,
	  	  tol_inner = 1e-8,
                  λ_approx=0,
		  v=randn(nep.n),
                  displaylevel=0,
		  max_it = 10,
		  eigsolver="default")


    levsolver = LinEigSolver();

	#The main program

	println("Running safeguarded iteration, initial approximation of λ: ",λ_approx)
	err=Inf;
	λ_m = λ_approx
	M = compute_Mder(nep,λ_m,0);
	v_m = v;
	d,v_m = levsolver.solve(compute_Mder(nep,λ_m,0),λ_t=λ_m,nev=1)
	residual = M*v_m;

	for k=1:max_it

		quotient_m = dot(v_m,(M*v_m));

		#Newton as an inner loop for the Rayleigh quotient
		println("k: ",k,"λ_m: ",λ_m);
		λ_m = compute_rf(nep,v_m,y=v_m, λ0=λ_m,TOL=tol_inner);
		println("k: ",k,"λ_m: ",λ_m);

		# Find closest eigenvalue to λ
        	# and the corresponding eigenvector
		#d,v_m = eigsolverfunc(nep,λ_m)
        d,v_m = levsolver.solve(compute_Mder(nep,λ_m,0),λ_t=λ_m,nev=1);

		# This to be changed ("errmeasure")		
		M = compute_Mder(nep,λ_m,0);
		residual = M*v_m;
		
		# Stopping criterion using err_measure()
		err=errmeasure(λ_m,v_m);
		if (err < tol_outer)
			println("Solution found at step ",k);
			println("λ = ",λ_m);
			return (λ_m,v_m)
		end

		println("Norm of residual at step ",k," = ",norm(residual));
		
		abs(d - λ_m)
		
	end
        		
 	err=errmeasure(λ_m,v_m)
	msg="Number of iterations exceeded. maxit=$(max_it)."
        throw(NoConvergenceException(λ_m,v_m,err,msg))

    end
    

end
