
# Newton-like methods for NEPs

    export newton
    export resinv
    export augnewton
    export quasinewton
    export newtonqr
    export implicitdet

#############################################################################
"""
    Newton raphsons method on nonlinear equation with (n+1) unknowns
"""
    newton(nep::NEP;params...)=newton(Complex128,nep;params...)
    function newton{T}(::Type{T},
                       nep::NEP;
                       errmeasure::Function =
                       default_errmeasure(nep::NEP),
                       tol=eps(real(T))*100,
                       maxit=10,
                       λ=zero(T),
                       v=randn(size(nep,1)),
                       c=v,
                       displaylevel=0,
                       armijo_factor=1,
                       armijo_max=5)

        # Ensure types λ and v are of type T
        λ=T(λ)
        v=Array{T,1}(v)
        c=Array{T,1}(c)

        err=Inf;
        v=v/dot(c,v);


        try
            for k=1:maxit
                err=errmeasure(λ,v)

                @ifd(print("Iteration:",k," errmeasure:",err))
                if (err< tol)
                    @ifd(print("\n"));
                    return (λ,v)
                end

                # Compute NEP matrix and derivative
                M = compute_Mder(nep,λ)
                Md = compute_Mder(nep,λ,1)

                # Create jacobian
                J = [M Md*v; c' 0];
                F = [M*v; c'*v-1];

                # Compute update
                delta=-J\F;  # Hardcoded backslash

                Δv=delta[1:size(nep,1)];
                Δλ=T(delta[size(nep,1)+1]);

                (Δλ,Δv,j,scaling)=armijo_rule(nep,errmeasure,err,
                                              λ,v,Δλ,Δv,real(T(armijo_factor)),armijo_max)
                if (j>0)
                    @ifd(@printf(" Armijo scaling=%f\n",scaling))
                else
                    @ifd(@printf("\n"))
                end


                # Update eigenvalue and eigvec
                v[:] += Δv
                λ = λ+Δλ
            end
        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))
            if (errmeasure(λ,v)>tol)
                # We need to compute an eigvec somehow
                v=compute_eigvec_from_eigval(nep,λ, default_linsolvercreator)
                v=v/dot(c,v)
            end
            return (λ,v)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end

    ############################################################################
"""
    Residual inverse iteration method for nonlinear eigenvalue problems.
"""
    resinv(nep::NEP;params...)=resinv(Complex128,nep;params...)
    function resinv{T}(::Type{T},
                       nep::NEP;
                       errmeasure::Function =
                       default_errmeasure(nep::NEP),
                       tol=eps(real(T))*100,
                       maxit=100,
                       λ=zero(T),
                       v=randn(real(T),size(nep,1)),
                       c=v,
                       displaylevel=0,
                       linsolvercreator::Function=default_linsolvercreator,
                       armijo_factor=1,
                       armijo_max=5)

        # Ensure types λ and v are of type T
        λ::T=T(λ)
        v=Array{T,1}(v)
        c=Array{T,1}(c)

        local linsolver::LinSolver=linsolvercreator(nep,λ)

        # If c is zero vector we take eigvec approx as left vector in
        # generalized Rayleigh functional
        use_v_as_rf_vector=false;
        if (norm(c)==0)
            use_v_as_rf_vector=true;
        end


        σ::T=λ;
        err=Inf;

        try
            for k=1:maxit
                # Normalize
                v = v/norm(v);

                err=errmeasure(λ,v)


                if (use_v_as_rf_vector)
                    c=v;
                end

                @ifd(@printf("Iteration: %2d errmeasure:%.18e ",k, err))
                @ifd(if (use_v_as_rf_vector); print(" v_as_rf_vector=",use_v_as_rf_vector); end)

                if (err< tol)
                    @ifd(print("\n"));
                    return (λ,v)
                end

                # Compute eigenvalue update
                λ1 = compute_rf(T, nep,v,y=c,λ0=λ,target=σ)
                Δλ=λ1-λ


                # Compute eigenvector update
                Δv = -lin_solve(linsolver,compute_Mlincomb(nep,λ1,v)) #M*v);


                (Δλ,Δv,j,scaling)=armijo_rule(nep,errmeasure,err,
                                              λ,v,Δλ,Δv,real(T(armijo_factor)),armijo_max)
                if (j>0 )
                    @ifd(@printf(" Armijo scaling=%f\n",scaling))
                else
                    @ifd(@printf("\n"))
                end

                # Update the eigenpair
                λ+=Δλ
                v[:] += Δv;

            end

        catch e

            if (!isa(e,Base.LinAlg.SingularException) && !isa(e,Base.LinAlg.LAPACKException))
                rethrow(e);
            end

            #isa(e, Base.LinAlg.SingularException) || ) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.

            @ifd(println("We have an exact eigenvalue."))
            if (errmeasure(λ,v)>tol)
                # We need to compute an eigvec somehow
                v= compute_eigvec_from_eigval(nep,λ, (nep, σ) -> linsolver)
                v=v/dot(c,v)
            end
            return (λ,v)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end





    # New augnewton
"""
    Augmented Newton's method. Equivalent to newton() but works only with operations on vectors of length n, instead of n+1.
"""
    augnewton(nep::NEP;params...)=augnewton(Complex128,nep;params...)
    function augnewton{T}(::Type{T},
                          nep::NEP;
                          errmeasure::Function = default_errmeasure(nep::NEP),
                          tol=eps(real(T))*100,
                          maxit=30,
                          λ=zero(T),
                          v=randn(real(T),size(nep,1)),
                          c=v,
                          displaylevel=0,
                          linsolvercreator::Function=backslash_linsolvercreator,
                          armijo_factor=1,
                          armijo_max=5)
        # Ensure types λ and v are of type T
        λ=T(λ)
        v=Array{T,1}(v)
        c=Array{T,1}(c)

        err=Inf;
        # If c is zero vector we take eigvec approx as normalization vector
        use_v_as_normalization_vector=false;
        if (norm(c)==0)
            use_v_as_normalization_vector=true;
            c = v /norm(v)^2
        end
        v=v/dot(c,v);
        local linsolver::LinSolver;

        try
            for k=1:maxit
                err=errmeasure(λ,v)
                @ifd(@printf("Iteration: %2d errmeasure:%.18e ",k, err))
                @ifd(if (use_v_as_normalization_vector); print(" v_as_normalization_vector=",use_v_as_normalization_vector); end)
                if (err< tol)
                    @ifd(print("\n"))
                    return (λ,v)
                end
                # tempvec =  (M(λ_k)^{-1})*M'(λ_k)*v_k
                # α = 1/(c'*(M(λ_k)^{-1})*M'(λ_k)*v_k);

                z=compute_Mlincomb(nep,λ,v,[T(1.0)],1)

                linsolver = linsolvercreator(nep,λ)
                tempvec = Array{T,1}(lin_solve(linsolver, z, tol=tol));

                if (use_v_as_normalization_vector)
                    c = v /norm(v)^2
                end
                α = T(1)/ dot(c,tempvec);

                Δλ=-α
                Δv=α*tempvec-v;

                (Δλ,Δv,j,scaling)=armijo_rule(nep,errmeasure,err,
                                              λ,v,Δλ,Δv,real(T(armijo_factor)),armijo_max)

                if (j>0)
                    @ifd(@printf(" Armijo scaling=%f\n",scaling))
                else
                    @ifd(@printf("\n"))
                end

                λ+=Δλ
                v+=Δv

            end

        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))
            if (errmeasure(λ,v)>tol)
                # We need to compute an eigvec
                v= compute_eigvec_from_eigval(nep,λ, linsolvercreator)
                v=v/dot(c,v)
            end
            return (λ,v)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end


"""
    quasinewton{T}([T=Complex128],nep::NEP,[errmeasure,][tol=eps(real(T))*100,][maxit=100,][λ=0,][v=randn(real(T),size(nep,1)),][ws=v,][displaylevel=0,][linsolvercreator::Function=default_linsolvercreator,][armijo_factor=1,][armijo_max=5])
An implementation of quasi-newton 2 as described in https://arxiv.org/pdf/1702.08492.pdf. The vector ws is a prepresentation of the normalization, in the sense that c'=ws'M(λ). The function has an implementation of armijo steplength control which may improve the convergence basin (or initial convergence), e.g., by setting `armijo_factor=0.5`.
"""
    quasinewton(nep::NEP;params...)=quasinewton(Complex128,nep;params...)
    function quasinewton{T}(::Type{T},
                           nep::NEP;
                           errmeasure::Function = default_errmeasure(nep::NEP),
                           tol=eps(real(T))*100,
                           maxit=100,
                           λ=zero(T),

                            v=randn(real(T),size(nep,1)),
                           ws=v,
                           displaylevel=0,
                           linsolvercreator::Function=default_linsolvercreator,
                           armijo_factor=1,
                           armijo_max=5)
        # Ensure types λ and v are of type T
        λ=T(λ)
        v=Array{T,1}(v)
        ws=Array{T,1}(ws) # Left vector such that c'=w'M(λ0) where c normalization

        err=Inf;

        local linsolver::LinSolver;
        @ifd(@printf("Precomputing linsolver (factorization)\n"))


        linsolver = linsolvercreator(nep,λ)

        try
            for k=1:maxit
                err=errmeasure(λ,v)
                @ifd(@printf("Iteration: %2d errmeasure:%.18e",k, err))
                @ifd(print(", λ=",λ))

                if (err< tol)
                    @ifd(@printf("\n"));
                    return (λ,v)
                end


                # Compute u=M(λ)v and w=M'(λ)v
                u=compute_Mlincomb(nep,λ,v,[T(1)],0);
                w=compute_Mlincomb(nep,λ,v,[T(1)],1);

                # Intermediate quantities
                Δλ=-dot(ws,u)/dot(ws,w);
                z=Δλ*w+u;
                Δv=-lin_solve(linsolver, z, tol=tol);


                (Δλ,Δv,j,scaling)=armijo_rule(nep,errmeasure,err,
                                              λ,v,Δλ,Δv,real(T(armijo_factor)),armijo_max)

                if (j>0)
                    @ifd(@printf(" Armijo scaling=%f\n",scaling));
                else
                    @ifd(@printf("\n"));
                end

                # Update eigenpair
                λ += Δλ
                v += Δv; # eigvec update

            end

        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))

            if (errmeasure(λ,v)>tol)
                # We need to compute an eigvec
                v= compute_eigvec_from_eigval(nep,λ, linsolvercreator)
                v=v/dot(c,v)
            end
            return (λ,v)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end


"""
    Newton-QR method.
"""
    newtonqr(nep::NEP;params...)=newtonqr(Complex128,nep;params...)
    function newtonqr{T}(::Type{T},
                       nep::NEP;
                       errmeasure::Function =
                       default_errmeasure(nep::NEP),
                       tol=eps(real(T))*100,
                       maxit=100,
                       λ=zero(T),
                       v=randn(real(T),size(nep,1)),
                       c=v,
                       displaylevel=0,
                       linsolvercreator::Function=default_linsolvercreator)


        n = size(nep,1);
        v = Array{T,1}(v);
        local err;

        en = zeros(n);
        en[n] = 1;
        try
            for k=1:maxit

                A = compute_Mder(nep,λ);
                Q,R,PI = qr(A,Val{true});#QR factorization with pivoting.

                P = eye(T,n)[:,PI];#The permutation matrix corresponding to the pivoted QR.

                p = R[1:n-1,1:n-1]\R[1:n-1,n];
                v = P*[-p;T(1)];#Right eigenvector
                w = Q*en;#Left eigenvector

                #err = abs(R[n,n])/vecnorm(compute_Mder(nep,λ),2);
                err=errmeasure(λ,v);
                @ifd(println("Iteration: ",k," errmeasure: ", err))
                if(err < tol)
                    return λ,v,w;
                end


                d = dot(Q[:,n],compute_Mlincomb(nep,λ,v,[T(1)],1));
                #d = dot(Q[:,n],compute_Mder(nep,λ,1)*P*[-p;T(1.0)]);
                λ = λ - R[n,n]/d;
            end
        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))

            if (k == maxit)
                # We need to compute an eigvec
                v= compute_eigvec_from_eigval(nep,λ, linsolvercreator)
                v=v/dot(c,v)
            end
            return (λ,v,w)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,err,msg))
    end


"""
    Implicit determinant method
"""
    implicitdet(nep::NEP;params...)=implicitdet(Complex128,nep;params...)
    function implicitdet{T}(::Type{T},
                       nep::NEP;
                       errmeasure::Function =
                       default_errmeasure(nep::NEP),
                       tol=eps(real(T))*100,
                       maxit=100,
                       λ=zero(T),
                       v=randn(real(T),size(nep,1)),
                       c=v,
                       displaylevel=0,
                       linsolvercreator::Function=default_linsolvercreator)


        n = size(nep,1);
        v = Array{T,1}(v);
        c = Array{T,1}(c);
        b = c;

        try
            for k=1:maxit

                A = compute_Mder(nep,λ);
                AA = [A b;c' 0];
                L,U,PI = lu(AA);


                P = eye(T,n+1)[PI,:];

                v = U\(L\(P*[zeros(T,n);T(1)]));
                #vp = U\(L\(P*[compute_Mlincomb(nep,λ,v[1:n],[T(-1.0)],1);0]));
                vp = U\(L\(P*[-1*compute_Mder(nep,λ,1)*v[1:n];0]))

                err = abs(v[n+1])/vecnorm(compute_Mder(nep,λ),2);
                @ifd(println("Iteration: ",k," errmeasure: ", err))
                if(err < tol)
                    @ifd(println(λ))
                    return λ,v[1:n];
                end

                λ = λ - v[n+1]/vp[n+1];
            end
        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            @ifd(println("We have an exact eigenvalue."))

            if (k == maxit)
                # We need to compute an eigvec
                v= compute_eigvec_from_eigval(nep,λ, linsolvercreator)
                v=v/dot(c,v)
            end
            return (λ,v)
        end

        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,v,NaN,msg))
    end

    # Armijo rule implementation
    function armijo_rule(nep,errmeasure,err0,λ,v,Δλ,Δv,armijo_factor,armijo_max)
        j=0
        if (armijo_factor<1)
            # take smaller and smaller steps until errmeasure is decreasing
            while (errmeasure(λ+Δλ,v+Δv)>err0 && j<armijo_max)
                j=j+1;
                Δv=Δv*armijo_factor;
                Δλ=Δλ*armijo_factor;
            end
        end
        return  (Δλ,Δv,j,armijo_factor^j)
end
