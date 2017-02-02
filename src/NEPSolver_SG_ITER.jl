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
                  λ_nr=1,
                  displaylevel=0,
		  max_it = 10)

	#The main loop

	println("Running safeguarded iteration, initial approximation of λ: ",λ_approx)

	λ_m = λ_approx

	M = compute_Mder(nep,λ_m,0);

	d,v_m = julia_eig(M);

	residual = M*v_m;
	quotient_m = dot(v_m,(M*v_m));


	for k=1:max_it

		quotient_m = dot(v_m,(M*v_m));

		#Newton as an inner loop (for the Rayleigh quotient)
		while (abs(quotient_m) > tol_inner)
			Md = compute_Mder(nep,λ_m,1);
			dλ = dot(v_m,(M*v_m))/dot(v_m,Md*v_m);
			println("dλ :",dλ)
			λ_m = λ_m - dλ;
			M = compute_Mder(nep,λ_m,0);
			quotient_m = dot(v_m,(M*v_m));
		end

		M = compute_Mder(nep,λ_m,0);
		# Find closest eigenvalue to λ
        	# and the corresponding eigenvector
		d,v_m = julia_eig(M);

		residual = M*v_m;

		if (norm(residual) < tol_outer)
			println("Solution found at step ",k);
			println("λ = ",λ_m);
			return (λ_m,v_m)
		end

		println("Norm of residual at step ",k," = ",norm(residual));
		
		abs(d - λ_m)
		
	end
        		
 	err=errmeasure(approx_m,approx_v)
	msg="Number of iterations exceeded. maxit=$(max_it)."
        throw(NoConvergenceException(approx_m,approx_v,err,msg))

    end



    #############################################################################
    #Call MATLAB eigs() 
    function matlab_eigs(A,λ = 0)


        aa=mxarray(A)
        s=mxarray(λ)

        @mput aa bb s
        @matlab begin
            s=double(s);
            aa=double(aa);
            (v,d)=eigs(aa,1,s);
        end
        @mget d v

        return d,v;
    end
    #############################################################################


    #Call Julia eig()
    function julia_eig(A,λ = 0)
        # Solve the linear eigenvalue problem
        D,V = eig(A);

        # Find closest eigenvalue to λ
        xx,idx=findmin(abs(D-λ))
        d=D[idx]

        # update eigenvector
        v=V[:,idx] 

        return d,v;     
    end

    

end
