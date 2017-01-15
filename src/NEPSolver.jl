module NEPSolver
    using NEPCore
    using MATLAB  # Dependence required for successive linear problems
    export newton
    export res_inv
    export mslp
    export aug_newton_old
    export aug_newton

    #############################################################################
    # Newton raphsons method on nonlinear equation with (n+1) unknowns
    function newton(nep::NEP;
                    errmeasure::Function =
                    default_errmeasure(nep::NEP),
                    tolerance=eps()*100,
                    maxit=10,
                    λ=0,
                    v=randn(nep.n),
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


    #############################################################################
    # Method of successive linear problems
    function mslp(nep::NEP;
                  errmeasure::Function =
                  default_errmeasure(nep::NEP),
                  tolerance=eps()*100,
                  maxit=100,
                  λ=0,
                  v=randn(nep.n),
                  c=v,
                  displaylevel=0,
                  eigsolver="default")

        σ=λ;     
        err=Inf;

        #Decide which solver will be called for the successive linear problems
        local eigsolverfunc::Function; 
        if(eigsolver == "eig")
            eigsolverfunc = julia_eig;

        elseif(eigsolver == "eigs") 
            println("Warning: Using eigsolver julia's eigs, which is not complete for generalized eigenvalue problems. See issue #1. ") # Issue #1
            eigsolverfunc = julia_eigs;

        elseif(eigsolver == "matlab_eigs")
            eigsolverfunc = matlab_eigs;

        else
            if issparse(compute_Mder(nep,λ,0)) 
                eigsolverfunc = matlab_eigs; # Default to matlab due to issue #1

            else
                eigsolverfunc = julia_eig;
            end
        end
            
            # Main loop
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

                # solve generalized eigenvalue problem
                d,v = eigsolverfunc(nep,λ,v);
                # update eigenvalue
                λ=λ-d

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
                v=(nep.Md(λ,0)+eps()*speye(nep.n))\v; # Requires matrix access
                v=v/dot(c,v)
            end
            return (λ,v)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
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
    function default_errmeasure(nep::NEP)
        f=function (λ,v);
            compute_resnorm(nep,λ,v)
        end
        return f
    end
        


    #############################################################################
    #Call MATLAB eigs() 
    function matlab_eigs(nep::NEP,λ = 0,v0=randn(nep.n))


        aa=mxarray(compute_Mder(nep,λ,0))
        bb=mxarray(compute_Mder(nep,λ,1))
        s=mxarray(λ)

        @mput aa bb s
        @matlab begin
            s=double(s);
            aa=double(aa);
            bb=double(bb);
            (v,d)=eigs(aa,bb,1,s);
        end
        @mget d v

        return d,v;
    end
    #############################################################################
    #Call Julia eig()
    function julia_eig(nep::NEP,λ = 0,v0=randn(nep.n))
        # Solve the linear eigenvalue problem
        D,V = eig(compute_Mder(nep,λ,0), compute_Mder(nep,λ,1));

        # Find closest eigenvalue to λ
        xx,idx=findmin(abs(D-λ))
        d=D[idx]

        # update eigenvector
        v=V[:,idx] 

        return d,v;     
    end

    #############################################################################
    #Call Julia eigs() 
    function julia_eigs(nep::NEP,λ = 0,v0=randn(nep.n))

        D,V = eigs(compute_Mder(nep,λ,0),nep,compute_Mder(nep,λ,1),
                   sigma=λ, v0=v0,nev=2,
                   tol=eps()*1000, maxiter=10)

        d=D[1]
        
        # update eigenvector
        v=V[:,1]

        return d,v;     
    end

end #End module
