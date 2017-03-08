



export jd_quad


function jd_quad(nep::NEP;
              errmeasure::Function =
              default_errmeasure(nep::NEP),
              tolerance=eps()*100,
              maxit=100,
              λ=0,
              v0=randn(nep.n),
              displaylevel=0,
              eigsolvertype::DataType=DefaultEigSolver)


   v = v0/norm(v0);
   
   V = zeros(nep.n,1);
   V[:,1] = v;
   u=v; 
   theta = λ; 

   #loop...

   for kk=1:maxit
    
    	err=errmeasure(theta,u)

    	if (displaylevel>0)
      		println("Iteration:",k," errmeasure:",err)
   	end
    	if (err< tolerance)
        	return (theta,u)
        end		

	# Projected matrices
	Ap0 = 	V'*(nep.A[1]*V);
	Ap1 = 	V'*(nep.A[2]*V);		
	Ap2 = 	V'*(nep.A[3]*V);	


	# Create the projected polynomial problem and 
	# find the eigenvalue with smallest absolute value

	Ap = [Ap0,Ap1,Ap2];
	pep_temp = PEP(Ap);
	Dc,Vc = polyeig(pep_temp,DefaultEigSolver);
	c = sortperm(abs(Dc));


	theta = Dc[c[1]];
	s = Vc[:,c[1]];
	s = s/norm(s);

	u = V*s;

	P = (eye(nep.n) - u*u');
	M = compute_Mder(nep,theta,0);
	r = M*u;

        X= (P*(M)*P);

	# Least squares, -pseudo_inv(X)*r
	Q,R = qr(X);	
	t = -R\(Q'*r);
  
    	#Modified Gram-Schmidt
    	for ii=1:kk
    	    temp = dot(V[:,ii],t);
    	    t = t - temp*V[:,ii];
    	end
    	# reorthogonalization 
    	for ii=1:kk
    	      temp = dot(V[:,ii],t);
    	      t = t - temp*V[:,ii];
    	end
    	v = t/norm(t);


	# Update the search space V
    	V = [V v];
    
    	
	println("Iteration: ",kk," norm of residual:", compute_resnorm(nep,theta,u))

    end

    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(λ,v,err,msg))
end

