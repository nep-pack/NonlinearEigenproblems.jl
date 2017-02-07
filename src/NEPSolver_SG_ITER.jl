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




	#Decide which solver will be called for the successive linear problems
        local eigsolverfunc::Function; 
        if(eigsolver == "eig")
            eigsolverfunc = julia_eig;
        elseif(eigsolver == "matlab_eigs")
            eigsolverfunc = matlab_eigs;
        else
            if issparse(compute_Mder(nep,λ,0)) 
                eigsolverfunc = matlab_eigs; # Default to matlab due to issue #1
            else
                eigsolverfunc = julia_eig;
            end
        end







	#The main program

	println("Running safeguarded iteration, initial approximation of λ: ",λ_approx)
	λ_m = λ_approx
	M = compute_Mder(nep,λ_m,0);
	v_m = v;
	d,v_m = eigsolverfunc(nep,λ_m);
	residual = M*v_m;

	for k=1:max_it

		quotient_m = dot(v_m,(M*v_m));

		#Newton as an inner loop for the Rayleigh quotient
		println("k: ",k,"λ_m: ",λ_m);
		λ_m = compute_rf(nep,v_m,y=v_m, λ0=λ_m,TOL=tol_inner);
		println("k: ",k,"λ_m: ",λ_m);

		#old inner loop
		#dλ=1.0;
		#while (abs(dλ) > tol_inner)
		#	Md = compute_Mder(nep,λ_m,1);
		#	dλ = dot(v_m,(M*v_m))/dot(v_m,Md*v_m);
		#	println("dλ :",dλ)
		#	λ_m = λ_m - dλ;
		#	M = compute_Mder(nep,λ_m,0);
		#	quotient_m = dot(v_m,(M*v_m));
		#end

		# Find closest eigenvalue to λ
        	# and the corresponding eigenvector
		d,v_m = eigsolverfunc(nep,λ_m)

		# This to be changed ("errmeasure")		
		M = compute_Mder(nep,λ_m,0);
		residual = M*v_m;

		if (norm(residual) < tol_outer)
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



    #############################################################################
    #Call MATLAB eigs() 
    function matlab_eigs(nep::NEP,λ = 0,v0=randn(nep.n))


        aa=mxarray(compute_Mder(nep,λ,0))
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
    function julia_eig(nep::NEP,λ = 0,v0=randn(nep.n))
        # Solve the linear eigenvalue problem
        D,V = eig(compute_Mder(nep,λ,0));

        # Find closest eigenvalue to λ
        xx,idx=findmin(abs(D-λ))
        d=D[idx]

        # update eigenvector
        v=V[:,idx] 

        return d,v;     
    end

    

end
