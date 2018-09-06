using LinearAlgebra
using Random

    export iar_old
    #Infinite Arnoldi for a given number of max iters (No error measure yet)
"""
    The Infinite Arnoldi method
"""
    function iar_old(
     nep::NEP;maxit=30,
     linsolvertype::DataType=DefaultLinSolver,tol=1e-12,Neig=maxit,
     errmeasure::Function = default_errmeasure(nep::NEP),
     σ=0.0,γ=1)

     n = nep.n; m = maxit;
     # initialization
     V = zeros(n*(m+1),m+1);
     H = zeros(m+1,m);
     y = zeros(n,m+1);
     α = [0;ones(m)];
     for i=2:m+1; α[i]=γ^(i-1); end	# rescaled coefficients
     local M0inv::LinSolver = linsolvertype(compute_Mder(nep,σ));
     err = zeros(m,m); # error history
     λ=complex(zeros(m+1)); Q=complex(zeros(n,m+1));

     V[1:n,1]=rand(n,1)/opnorm(randn(n,1));

     k=1; conv_eig=0;
     while (k <= m)&(conv_eig<=Neig)

      #y[:,2:k+1] = reshape(V[1:n*k,k],n,k);

      y[:,2:k+1] = reshape(view(V,1:1:n*k,k),n,k);
      # no improvement, just sperimenting.
      for j=1:k
       y[:,j+1]=y[:,j+1]/j;
      end

      y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],a=α[1:k+1]);
      y[:,1] = -lin_solve(M0inv,y[:,1]);

      vv=reshape(y[:,1:k+1],(k+1)*n,1);
      # orthogonalization
      h,vv = doubleGS(V,vv,k,n);
      H[1:k,k]=h;
      beta=opnorm(vv);

      H[k+1,k]=beta;
      V[1:(k+1)*n,k+1]=vv/beta;

      # compute error history
      D,Z = eigen(H[1:k,1:k]); D=σ+γ./D;

      conv_eig=0;
      for s=1:k
       err[k,s]=errmeasure(D[s],V[1:n,1:k]*Z[:,s]);
       if err[k,s]>10; err[k,s]=1; end	# artificial fix
       if err[k,s]<tol
        conv_eig=conv_eig+1;
        Q[:,conv_eig]=V[1:n,1:k]*Z[:,s]; λ[conv_eig]=D[s];
       end
      end

      k=k+1;
      end

      # extract the converged Ritzpairs
      λ=λ[1:min(length(λ),conv_eig)];
      Q=Q[:,1:min(size(Q,2),conv_eig)];

      return λ,Q,err
    end


    function doubleGS(V,vv,k,n)

            h=V[1:(k+1)*n,1:k]'*vv;

            vv=vv-V[1:(k+1)*n,1:k]*h;

            g=V[1:(k+1)*n,1:k]'*vv;
            vv=vv-V[1:(k+1)*n,1:k]*g;

            h = h+g;
            return h,vv;
    end

    function singleGS(V,vv,k,n)

            h=V[1:(k+1)*n,1:k]'*vv;

            vv=vv-V[1:(k+1)*n,1:k]*h;

            return h,vv;
    end
