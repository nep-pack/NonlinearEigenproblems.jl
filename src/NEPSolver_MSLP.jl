"""
    Module with Method of successive linear problems
"""
module NEPSolver_MSLP

    using NEPCore
    using MATLAB  # Dependence required for successive linear problems
    export mslp



    #############################################################################
    # Method of successive linear problems
    function mslp(nep::NEP;
                  errmeasure::Function =
                  default_errmeasure(nep::NEP),
                  tolerance=eps()*100,
                  maxit=100,
                  λ=0,
                  v=randn(nep.n),
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
        elseif(eigsolver == "inv_it")
            eigsolverfunc = inv_it
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
            v=v/norm(v);

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

        D,V = eigs(compute_Mder(nep,λ,0),compute_Mder(nep,λ,1),
                   sigma=0,nev=2,
                   tol=eps()*1000, maxiter=10)

        d=D[1]

        # update eigenvector
        v=V[:,1]

        return d,v;     
    end

    # A naive implementation of inverse iteration of generalized
    # linear eigenvalue problems
    function inv_it(nep::NEP,λ=0,v0=randn(nep.n),iters=2)
        Mp=compute_Mder(nep,λ,1);
        M=compute_Mder(nep,λ)
        local w=copy(v0)
        for i=1:iters
            w=M\(Mp*w)
            w=w/norm(w)
        end
        Δ=dot(w,M*w)/dot(w,Mp*w) # Comp delta with Rayleigh Quotient
        return Δ,w
    end

end
