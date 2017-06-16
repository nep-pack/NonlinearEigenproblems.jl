export rfi
"""
Two-sided Rayleigh functional Iteration, as given as Algorithm 4 in  "Nonlinear Eigenvalue Problems: Newton-type Methods and
Nonlinear Rayleigh Functionals", by Kathrin Schreiber.
"""
function rfi(nep::NEP;
            errmeasure::Function=default_errmeasure(nep::NEP),
            tolerance = eps()*100,
            maxit=100,
            λ = 0.0,
            v = randn(nep.n),
            u = randn(nep.n),
            linsolvercreator::Function=default_linsolvercreator,
            displaylevel=1)

        err = Inf;

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

                #S1: T(λ_k)x_(k+1) = T'(λ_k)u_(k)
                local linsolver::LinSolver = linsolvercreator(nep,λ);
                x = lin_solve(linsolver,compute_Mder(nep,λ,1)*u,tol = tolerance);
                u = x/norm(x);

                #S2: (T(λ_k)^H)y_(k+1) = (T'(λ_k)^H)v_(k)
                L = compute_Mder(nep,λ)';
                R = compute_Mder(nep,λ,1)'*v;
                y = L\R;#Hardcoded backslash
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