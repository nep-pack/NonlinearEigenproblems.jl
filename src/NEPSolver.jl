module NEPSolver
    using NEPCore
    export newton
    export res_inv
    export aug_newton_old
    export aug_newton
    export iar

    #############################################################################
    # Newton raphsons method on nonlinear equation with (n+1) unknowns
    function newton(nep::NEP;
                    errmeasure::Function =
                    default_errmeasure(nep::NEP),
                    tolerance=eps()*100,
                    maxit=10,
                    λ=0,
                    v=randn(size(nep,1)),
                    c=v,
                    displaylevel=0)

        err=Inf;
        v=v/dot(c,v);

        try
            for k=1:maxit
                err=errmeasure(λ,v)
                if (displaylevel>0)
                    println("Iteration:",k," errmeasure:",err)
                end
                if (err< tolerance)
                    return (λ,v)
                end

                # Compute NEP matrix and derivative
                M=compute_Mder(nep,λ)
                Md=compute_Mder(nep,λ,1)

                # Create jacobian
                J=[M Md*v; c' 0];
                F=[M*v; c'*v-1];

                # Compute update
                delta=-J\F;

                # Update eigenvalue and eigvec
                v=v+delta[1:nep.n];
                λ=λ+delta[nep.n+1];

            end
        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            if (displaylevel>0)
                println("We have an exact eigenvalue.")
            end
            if (errmeasure(λ,v)>tolerance)
                # We need to compute an eigvec somehow
                v=(compute_Mder(nep,λ,0)+eps()*speye(nep.n))\v; # Requires matrix access
                v=v/dot(c,v)
            end
            return (λ,v)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end


    ############################################################################
    function res_inv(nep::NEP;
                     errmeasure::Function =
                     default_errmeasure(nep::NEP),
                     tolerance=eps()*100,
                     maxit=100,
                     λ=0,
                     v=randn(nep.n),
                     c=v,
                     displaylevel=0,
                     linsolver=LinSolver(compute_Mder(nep,λ))
                     )

        σ=λ;
        # Compute a (julia-selected) factorization of M(σ)
        #Mσ=factorize(nep.Md(σ));
        #linsolver=LinSolver(nep.Md(σ))

        err=Inf;
        try
            for k=1:maxit
                # Normalize
                v=v/dot(c,v);


                err=errmeasure(λ,v)


                if (displaylevel>0)
                    println("Iteration:",k," errmeasure:",err)
                end
                if (err< tolerance)
                    return (λ,v)
                end

                # Compute eigenvalue update
                λ=compute_rf(nep,v,y=c,λ0=λ,target=σ)

                # Re-compute NEP matrix and derivative
                M=compute_Mder(nep,λ);

                # Compute eigenvector update
	        # Δv=Mσ\nep.Mlincomb(λ,v) #M*v);
                tol=eps()
	        Δv=linsolver.solve(compute_Mlincomb(nep,λ,reshape(v,nep.n,1)),
                                   tol=tol) #M*v);

                # Update the eigvector
                v=v-Δv;

            end

        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            if (displaylevel>0)
                println("We have an exact eigenvalue.")
            end
            if (errmeasure(λ,v)>tolerance)
                # We need to compute an eigvec somehow
                v=(nep.Md(λ,0)+eps()*speye(nep.n))\v; # Requires matrix access
                v=v/dot(c,v)
            end
            return (λ,v)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end


       


    # New aug_newton 
    function aug_newton(nep::NEP;
                        errmeasure::Function = default_errmeasure(nep::NEP),
                        tolerance=eps()*100,
                        maxit=30,
                        λ=0,
                        v=randn(nep.n),
                        c=v,
                        displaylevel=0)
        
        err=Inf;
        v=v/dot(c,v);
        try
            for k=1:maxit
                #err=errmeasure(λ,v)
                err=errmeasure(λ,reshape(v,nep.n))
                if (displaylevel>0)
                    println("Iteration:",k," errmeasure:",err)
                end
                if (err< tolerance)
                    return (λ,v)
                end
                # tempvec =  (M(λ_k)^{-1})*M'(λ_k)*v_k
                # α = 1/(c'*(M(λ_k)^{-1})*M'(λ_k)*v_k);
                z=compute_Mlincomb(nep,λ,v*ones(1,2),a=[0,1])

                tempvec = compute_Mder(nep,λ)\z
                α = 1.0/dot(c,tempvec);

                λ = λ - α;

                v = α*tempvec;

            end

        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            if (displaylevel>0)
                println("We have an exact eigenvalue.")
            end
            if (errmeasure(λ,v)>tolerance)
                # We need to compute an eigvec somehow
                v=(compute_Mder(nep,λ,0)+eps()*speye(nep.n))\v; # Requires matrix access
                v=v/dot(c,v)
            end
            return (λ,v)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end
      
    #Infinite Arnoldi for a given number of max iterations(No error measure yet)  
    function iar(nep::NEP,maxit=30)

        n = nep.n;            

        m = maxit;

        #####Pre-allocating the necessary matrices and vectors######
        V = zeros(n*(m+1),m+1);#Krylov subspace

        H = zeros(m+1,m);#Hessenberg matrix

        Bv = zeros(n,m+1);#For storing y0,y1,y2,.....,y_{k+1} at the kth iteration,
                        #where V_{k+1} = vec(y0,y1,y2,....,y_{k+1})

        α = [0;ones(m)];#Coefficients for the LC: 0*M(0)*y0+∑M^{i}(0)*y_i

        W = zeros(n,m+1);#For storing W[:,1:k+1] = (0,y1,2*y2,3*y3,.......,k*y_k)

        M0inv = LinSolver(compute_Mder(nep,0.0));#For computing the action of M(0)^{-1} later by M0inv.solve()

        V[1:n,1]=rand(n,1)/norm(randn(n,1));#Initializing the basis

        for k=1:m
            ########## Compute action of the operator B in Bv #########

            #Compute y0 = Bv[1:n,1]  
            W[:,2:k+1] = reshape(V[1:n*k,k],n,k);#Extract v_{k+1} and reshape it into a matrix   
            Bv[1:n,1] = compute_Mlincomb_from_Mder(nep,0.0,W,α[1:k+1]);
            Bv[1:n,1] = -M0inv.solve(Bv[1:n,1]);

            #Compute y1,y2,......y_k
            for j=2:k+1
                Bv[:,j]=W[:,j]/j; #Vectorizable 
            end
 
            #vv = V_{k+1} = vec(y0,y1,y2,....,y_{k+1}
            vv=reshape(Bv[:,1:k+1],(k+1)*n,1);

            # double GS-orthogonalization
            h,vv = doubleGS(V,vv,k,n);
            H[1:k,k]=h;

            beta=norm(vv);

            H[k+1,k]=beta;

            V[1:(k+1)*n,k+1]=vv/beta;
        end


        D,V=eig(H[1:m,1:m]);
        
        D=1./D;

        return D,V
    end

    function doubleGS(V,vv,k,n)

            h=V[1:(k+1)*n,1:k]'*vv;

            vv=vv-V[1:(k+1)*n,1:k]*h;
 
            g=V[1:(k+1)*n,1:k]'*vv;
            vv=vv-V[1:(k+1)*n,1:k]*g;

            h = h+g;
            return h,vv;
    end
    ##############################################################################
    #  function default_errmeasure(nep::NEP, displaylevel)
    #      # If no relresnorm available use resnorm
    #      if (isdefined(nep, :relresnorm))
    #          return nep.relresnorm;
    #      else
    #          if (displaylevel>0)
    #              println("Using resnorm")
    #          end
    #          return nep.resnorm;
    #      end
    #  end
    #
        



end #End module
