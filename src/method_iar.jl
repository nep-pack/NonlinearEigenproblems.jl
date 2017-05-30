    export iar
    #Infinite Arnoldi for a given number of max iters (No error measure yet) 
"""
    The Infinite Arnoldi method 
"""
    function iar(
     nep::NEP;
     maxit=30,
     linsolvercreator::Function=default_linsolvercreator,
     tol=1e-12,
     Neig=maxit,                                  
     errmeasure::Function = default_errmeasure(nep::NEP),
     σ=0.0,
     γ=1,
     displaylevel=0
    )

     n = size(nep,1); m = maxit;
     # initialization
     V = zeros(Complex128,n*(m+1),m+1);
     H = zeros(Complex128,m+1,m);
     y = zeros(Complex128,n,m+1);
     α = [0;ones(m)];
     # rescaled coefficients(TODO: integrate in compute_Mlincomb)
     for i=2:m+1; α[i]=γ^(i-1); end 
     local M0inv::LinSolver = linsolvercreator(nep,σ);
     err = zeros(m,m); 			
     λ=complex(zeros(m+1)); Q=complex(zeros(n,m+1));

     vv=view(V,1:1:n,1); # next vector V[:,k+1]
     vv[:]=rand(n,1); vv[:]=vv[:]/norm(vv);

     k=1; conv_eig=0;
        while (k <= m)&(conv_eig<=Neig)
            if (displaylevel>0)
                println("Iteration:",k)
                end
      VV=view(V,1:1:n*(k+1),1:k); # extact subarrays, memory-CPU efficient
      vv=view(V,1:1:n*(k+1),k+1); # next vector V[:,k+1]

      y[:,2:k+1] = reshape(VV[1:1:n*k,k],n,k);
      for j=1:k	
       y[:,j+1]=y[:,j+1]/j;  
      end

      y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],a=α[1:k+1]);
      y[:,1] = -lin_solve(M0inv,y[:,1]);

      vv[:]=reshape(y[:,1:k+1],(k+1)*n,1);
      # orthogonalization
      h,vv[:] = doubleGS(VV,vv,k,n);
      H[1:k,k]=h;
      beta=norm(vv);

      H[k+1,k]=beta;
      vv[:]=vv[:]/beta;

      # compute error history
      D,Z=eig(H[1:k,1:k]); D=σ+γ./D;
      VV=view(V,1:1:n,1:k);	# extract proper subarray 
      conv_eig=0;
      for s=1:k
       err[k,s]=errmeasure(D[s],VV*Z[:,s]);
       if err[k,s]>10; err[k,s]=1; end	# artificial fix
       if err[k,s]<tol
        conv_eig=conv_eig+1;
        Q[:,conv_eig]=VV*Z[:,s]; λ[conv_eig]=D[s];
       end
      end

      k=k+1;
      end

      # extract the converged Ritzpairs
      λ=λ[1:min(length(λ),conv_eig)];
      Q=Q[:,1:min(size(Q,2),conv_eig)];

      return λ,Q,err
    end


    function doubleGS(VV,vv,k,n)

            h=VV'*vv;

            vv=vv-VV*h;
 
            g=VV'*vv;
            vv=vv-VV*g;

            h = h+g;
            return h,vv;
    end

    function singleGS(V,vv,k,n)

            h=V[1:(k+1)*n,1:k]'*vv;

            vv=vv-V[1:(k+1)*n,1:k]*h;

            return h,vv;
    end
