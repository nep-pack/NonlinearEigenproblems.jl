using IterativeSolvers
using LinearAlgebra
using Random

export nlar
export default_eigval_sorter
export residual_eigval_sorter
export threshold_eigval_sorter



"""
    function nlar([eltype],nep::ProjectableNEP,[orthmethod=ModifiedGramSchmidt()],[neigs=10],[errmeasure],[tol=eps(real(T))*100],[maxit=100],[0=0],[v=randn(T,size(nep,1))],[logger=0],[linsolvercreator=DefaultLinSolverCreator()],[R=0.01],[eigval_sorter=residual_eigval_sorter],[qrfact_orth=false],[max_subspace=100],[num_restart_ritz_vecs=8],[inner_solver_method=DefaultInnerSolver(),][inner_logger=0])

The function implements the Nonlinear Arnoldi method, which finds `neigs` eigenpairs (or throws a `NoConvergenceException`) by projecting the problem to a subspace that is expanded in the course  of the algorithm.
The basis is orthogonalized either by using the QR method if `qrfact_orth` is `true` or else by an orthogonalization method `orthmethod`).
This entails solving a smaller projected problem using a method specified by `inner_solver_method`.
The logging of the inner solvers are descided by `inner_logger`, which works in the same way as `logger`.
(`λ`,`v`) is the initial guess for the eigenpair. `linsolvercreator` specifies how the linear system is created and solved.
`R` is a parameter used by the function specified by `eigval_sorter` to reject those ritz values that are within a distance `R` from any of the converged eigenvalues, so that repeated convergence to the same eigenpair can be avoided.
`max_subspace` is the maximum allowable size of the basis befor the algorithm restarts using a basis made of `num_restart_ritz_vecs` ritz vectors and the eigenvectors that the algorithm has converged to.

See [`augnewton`](@ref) for other parameters.


# Example
```julia-repl
julia> nep=nep_gallery("dep0_tridiag");
julia> λ,v=nlar(nep,tol=1e-7,neigs=1,maxit=100,v=ones(size(nep,1)));
julia> norm(compute_Mlincomb(nep,λ[1],v))
8.00192341259751e-7
```

# References
* H. Voss, An Arnoldi method for nonlinear eigenvalue problems. BIT. Numer. Math. 44: 387-401, 2004.

"""
nlar(nep::NEP;params...) = nlar(ComplexF64,nep::NEP;params...)
function nlar(::Type{T},
            nep::ProjectableNEP;
            orthmethod = ModifiedGramSchmidt(),
            neigs::Int=10,                                     #Number of eigenvalues required
            errmeasure::ErrmeasureType = DefaultErrmeasure(nep),
            tol = eps(real(T))*100,
            maxit::Int = 100,
            λ::Number = zero(T),
            v::Vector = randn(T,size(nep,1)),
            logger = 0,
            linsolvercreator=DefaultLinSolverCreator(),
            R = 0.01,
            eigval_sorter::Function = residual_eigval_sorter, #Function to sort eigenvalues of the projected NEP
            qrfact_orth::Bool = false,
            max_subspace::Int = 100,                           #Maximum subspace size before we implement restarting
            num_restart_ritz_vecs::Int=8,
            inner_solver_method = DefaultInnerSolver(),
            inner_logger = 0) where {T<:Number}

        @parse_logger_param!(logger)
        @parse_logger_param!(inner_logger)

        #Check if maxit is larger than problem size
        if (maxit > size(nep,1))
            @warn "Maximum iteration count maxit=$maxit larger than problem size n=$(size(nep,1)). Reducing maxit."
            maxit = size(nep,1);
        end

        #Check if number of ritz vectors for restating is greater than neigs
        if(num_restart_ritz_vecs > neigs)
            @warn "Number of ritz vectors for restarting num_restart_ritz_vecs=$num_restart_ritz_vecs larger than neigs=$neigs. Reducing num_restart_ritz_vecs."
            num_restart_ritz_vecs = neigs;
        end

        #Check if maximum allowable subspace size is lesser than num_restart_ritz_vecs
        if(max_subspace < num_restart_ritz_vecs)
            @warn "Maximum subspace max_subspace=$max_subspace smaller than number of restarting ritz vectors num_restart_ritz_vecs=$num_restart_ritz_vecs. Increasing max_subspace"
            max_subspace = num_restart_ritz_vecs+20; #20 is hardcoded.
        end

        local σ::T = T(λ); #Initial pole
        n = size(nep,1)
        λ::T = T(λ)
        nu::T = λ
        u::Vector{T} = v

        #Initialize the basis V_1
        V::Matrix{T}= zeros(T,n,max_subspace);
        X::Matrix{T} = zeros(T,n,neigs);
        V[:,1] = normalize(v);
        cbs = 1;#Current basis size


        D::Vector{T} = zeros(T,neigs);#To store the converged eigenvalues
        err_hist=eps()*ones(maxit,neigs) ;#Error history

        Z::Matrix{T} = zeros(T,n,neigs+num_restart_ritz_vecs);#The matrix used for constructing the restarted basis
        m = 0; #Number of converged eigenvalues

        k = 1;

        proj_nep = create_proj_NEP(nep,maxit,T);

        local linsolver::LinSolver=create_linsolver(linsolvercreator,nep,λ)

        err = Inf;



        push_info!(logger, "Using inner solver $inner_solver_method");

        while ((m < neigs) && (k < maxit))
            Vk = view(V,:,1:cbs)
            # Construct and solve the small projected PEP projected problem (V^H)T(λ)Vx = 0
            expand_projectmatrices!(proj_nep,Vk,Vk);

            #Use inner_solve() to solve the smaller projected problem
            push_info!(logger, 2, "Solving inner problem",continues=true)
            dd,vv = inner_solve(inner_solver_method,T,proj_nep,neigs=neigs,σ=σ,inner_logger=inner_logger);
            push_info!(logger, 2, ". Done.")
            # Sort the eigenvalues of the projected problem
            nuv,yv = eigval_sorter(nep,dd,vv,σ,D,R,Vk);
            nu = nuv[1]; y=yv[:,1];

            if (isinf(nu))
                error("We did not find any (non-converged) eigenvalues to target")
            end


            #Determine ritz vector and residual
            u[:] = Vk*y; # Note: y and u are vectors (not matrices)

            #Normalize and compute residual
            normalize!(u);
            res::AbstractVector = compute_Mlincomb(nep,nu,u);


            #Check for convergence of one of the eigenvalues
            err = estimate_error(errmeasure,nu,u);


            push_iteration_info!(logger,k,λ=nu,v=u,err=err);
            err_hist[k,m+1]=err;
            if(err < tol)
                mplusone=m+1;
                push_info!(logger,
                           "****** $mplusone converged to eigenvalue: $nu errmeasure:$err");

                #Add to the set of converged eigenvalues and eigenvectors
                D[m+1] = nu;
                X[:,m+1] = u;

                ## Sort and select he eigenvalues of the projected problem as described before
                nuv,yv = eigval_sorter(nep,dd,vv,σ,D,R,Vk);
                nu1=nuv[1];
                y1=yv[:,1];

                #Compute residual again
                u1 = Vk*y1;
                normalize!(u1);
                res = compute_Mlincomb(nep,nu1,u1);

                m = m+1;
            end

            #Check if basis size has exceeded max_subspace. If yes, then restart.
            if(size(Vk,2) >= max_subspace)
                cbs = m+num_restart_ritz_vecs;

                #Construct the new basis
                Zv = view(Z,:,1:cbs);
                Zv[:,1:m] = X[:,1:m]; #Set the first m vectors of the restarted basis to be the converged eigenvectors
                Zv[:,m+1:cbs] = Vk*yv[:,1:num_restart_ritz_vecs]; #Set the rest to be the ritz vectors of the projected problem

                Q,_ = qr(Zv); #Thin QR-factorization
                V[:,1:cbs] = Matrix(Q);

            else
                #Compute new vector Δv to add to the search space V(k+1) = (Vk,Δv)
                Δv=lin_solve(linsolver,res)
                #Orthogonalize and normalize
                if (qrfact_orth)
                    # Orthogonalize the entire basis matrix
                    # together with Δv using QR-method.
                    # Slow but robust.
                    Q,_ = qr(hcat(Vk,Δv))
                    Q = Matrix(Q)

                    cbs = cbs+1;
                    V[:,1:cbs]=Q;
                else
                    h=zeros(T,k);
                    orthogonalize_and_normalize!(Vk,Δv,h,orthmethod);

                    #Expand basis
                    cbs = cbs+1;
                    V[:,cbs] = Δv;
                end
            end
            k = k+1;
        end


        #Throw no convergence exception if enough eigenvalues were not found.
        if(k >= maxit && m < neigs)
            msg="Number of iterations exceeded. maxit=$(maxit) and only $(m) eigenvalues converged out of $(neigs)."
            throw(NoConvergenceException(nu,u,err,msg))
        end

        return D,X,err_hist;
    end


###############################################################################################################
# Ritz value discarder:
# This function is called as a first step by all the sorter functions to discard ritz values that lie within
# a certain distance R of the converged eigenvalues.
function discard_ritz_values!(dd,D,R)
    for i=1:size(dd,1)
        for j=1:size(D,1)
            if (abs(dd[i]-D[j])<R)
                dd[i]=Inf; #Discard all Ritz values within a particular radius R
            end
        end
    end

end



# Default ritzvalue sorter:
# First discard all Ritz values within a distance R of any of the converged eigenvalues(of the original problem).
# Then sort by distance from the shift and select the mm-th furthest value from the pole.

## D = already computed eigenvalues
## dd, vv eigenpairs of projected problem
## σ targets

function  default_eigval_sorter(nep::NEP,dd,vv,σ,D,R,Vk)
    dd2=copy(dd);

    #Discard ritz values within a distance R of the converged eigenvalues
    discard_ritz_values!(dd2,D,R)

    ii = sortperm(abs.(dd2-σ));

    nu = dd2[ii];
    y = vv[:,ii];

    return nu,y;
end


# Residual-based Ritz value sorter:
# First discard all Ritz values within a distance R of any of the converged eigenvalues(of the original problem).
# Then select that Ritz value which gives the mm-th minimum product of (residual and distance from pole).
function residual_eigval_sorter(nep::NEP,dd,vv,σ,D,R,Vk,errmeasure::ErrmeasureType = DefaultErrmeasure(nep))

    eig_res = zeros(size(dd,1));
    dd2=copy(dd);

    #Discard ritz values within a distance R of the converged eigenvalues
    discard_ritz_values!(dd2,D,R)


    #Compute residuals for each Ritz value
    for i=1:size(dd,1)
        eig_res[i] = estimate_error(errmeasure,dd[i],Vk*vv[:,i]);
    end

    #Sort according to methods
    ii = sortperm(eig_res .* abs.(dd2 .- σ))

    nu = dd[ii];
    y = vv[:,ii];

    return nu,y;
end

# Threshold residual based Ritz value sorter:
# Same as residual_eigval_sorter() except that errors above a certain threshold are set to the threshold.
function threshold_eigval_sorter(nep::NEP,dd,vv,σ,D,R,Vk,errmeasure::ErrmeasureType = DefaultErrmeasure(nep),threshold=0.1)

    eig_res = zeros(size(dd,1));
    dd2=copy(dd);

    #Discard ritz values within a distance R of the converged eigenvalues
    discard_ritz_values!(dd2,D,R)

    #Compute residuals for each Ritz value
    temp_res = 0;
    for i=1:size(dd,1)
        temp_res = error_measure(errmeasure,dd[i],Vk*vv[:,i]);
        if(temp_res > threshold)
            eig_res[i] = threshold;
        else
            eig_res[i] = temp_res;
        end
    end

    #Sort according to methods
    ii = sortperm(eig_res.*abs.(dd2-σ));

    nu = dd[ii];
    y = vv[:,ii];

    return nu,y;
end
