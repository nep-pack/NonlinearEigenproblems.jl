export rfi
export rfi_b
"""
Two-sided Rayleigh functional Iteration, as given as Algorithm 4 in  "Nonlinear Eigenvalue Problems: Newton-type Methods and
Nonlinear Rayleigh Functionals", by Kathrin Schreiber.
"""
function rfi(nep::NEP,
            nept::NEP;
            errmeasure::Function=default_errmeasure(nep::NEP),
            tolerance = eps()*1000,
            maxit=100,
            λ = 0.0+0.0im,
            v = randn(nep.n),
            u = randn(nep.n),
            linsolvercreator::Function=default_linsolvercreator,
            displaylevel=1)

        err = Inf;

        #Ensure type coherence
        T = typeof(λ);
        v = Array{T,1}(v);
        u = Array{T,1}(u);
        #Normalize v and u
        v = v/norm(v);
        u = u/norm(u);


        try
            for k=1:maxit
                err = errmeasure(λ,u);

                if(err < tolerance)
                    return λ,u
                end

                if (displaylevel>0)
                    @printf("Iteration: %2d errmeasure:%.18e \n",k, err);
                end

                local linsolver::LinSolver = linsolvercreator(nep,λ);
                local linsolver_t::LinSolver = linsolvercreator(nept,λ);
                
                #S1: T(λ_k)x_(k+1) = T'(λ_k)u_(k) 
                x = lin_solve(linsolver,compute_Mlincomb(nep,λ,u,[1],1),tol = tolerance);
                u = x/norm(x);

                #S2: (T(λ_k)^H)y_(k+1) = (T'(λ_k)^H)v_(k)
                y = lin_solve(linsolver_t,compute_Mlincomb(nept,λ,v,[1],1),tol = tolerance);
                v = y/norm(y);

                λ =compute_rf(nep,u;y=v);
            end
        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            if (displaylevel>0)
                println("We have an exact eigenvalue.")
            end
            if (errmeasure(λ,u)>tolerance)
                # We need to compute an eigvec somehow
                u=(nep.Md(λ,0)+eps(real(T))*speye(size(nep,1)))\u; # Requires matrix access
                u=u/norm(u);
            end
            return (λ,u)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,u,err,msg))
end

"""
Two-sided Rayleigh functional Iteration(Bordered version), as given as Algorithm 5 in  "Nonlinear Eigenvalue Problems: Newton-type Methods and
Nonlinear Rayleigh Functionals", by Kathrin Schreiber.
"""
function rfi_b(nep::NEP,
            nept::NEP;
            errmeasure::Function=default_errmeasure(nep::NEP),
            tolerance = eps()*1000,
            maxit=100,
            λ = 0.0+0.0im,
            v = randn(nep.n),
            u = randn(nep.n),
            linsolvercreator::Function=default_linsolvercreator,
            displaylevel=1)

        err = Inf;

        #Ensure type coherence
        T = typeof(λ);
        v = Array{T,1}(v);
        u = Array{T,1}(u);
        #Normalize v and u
        v = v/norm(v);
        u = u/norm(u);


        try
            for k=1:maxit
                err = errmeasure(λ,u);

                if(err < tolerance)
                    return λ,u
                end

                if (displaylevel>0)
                    @printf("Iteration: %2d errmeasure:%.18e \n",k, err);
                end

                #Construct C_k
                C = [compute_Mder(nep,λ,0) compute_Mlincomb(nep,λ,u,[1],1);v'*compute_Mder(nep,λ,1) T(0.0)];
                #local linsolver::LinSolver = BackslashLinSolver(C);

                #C[s;μ] = -[T(λ)u;0]
                #l1 = lin_solve(linsolver,-[compute_Mlincomb(nep,λ,u,[1],0);0],tol = tolerance);
                l1 = C\-[compute_Mlincomb(nep,λ,u,[1],0);0];
                s = l1[1:end-1];
                u = (u+s)/norm(u+s);
                
                #C[t;ν] = -[T(λ)'v;0]
                #l2 = lin_solve(linsolver,-[compute_Mlincomb(nept,λ,v,[1],0);0],tol = tolerance);
                l2 = C\-[compute_Mlincomb(nept,λ,v,[1],0);0];
                t = l2[1:end-1];
                v = (v+t)/norm(v+t);

                λ =compute_rf(nep,u;y=v);
            end
        catch e
            isa(e, Base.LinAlg.SingularException) || rethrow(e)
            # This should not cast an error since it means that λ is
            # already an eigenvalue.
            if (displaylevel>0)
                println("We have an exact eigenvalue.")
            end
            if (errmeasure(λ,u)>tolerance)
                # We need to compute an eigvec somehow
                u=(nep.Md(λ,0)+eps(real(T))*speye(size(nep,1)))\u; # Requires matrix access
                u=u/norm(u);
            end
            return (λ,u)
        end
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,u,err,msg))
end