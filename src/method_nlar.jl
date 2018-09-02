using IterativeSolvers
using LinearAlgebra
using Random

export nlar
export default_eigval_sorter
export residual_eigval_sorter
export threshold_eigval_sorter

###############################################################################################################
# Default ritzvalue sorter:
# First discard all Ritz values within a distance R of any of the converged eigenvalues(of the original problem).
# Then sort by distance from the shift and select the mm-th furthest value from the pole.

## D = already computed eigenvalues
## dd, vv eigenpairs of projected problem
## σ targets
"""
 The Nonlinear Arnoldi method, as introduced in "An Arnoldi method for nonlinear eigenvalue problems" by H.Voss
"""
function  default_eigval_sorter(nep::NEP,dd,vv,σ,D,mm,R,Vk)
    dd2=copy(dd);

    ## Check distance of each eigenvalue of the projected NEP(i.e. in dd)
    ## from each eigenvalue that as already converged(i.e. in D)
    for i=1:size(dd,1)
        for j=1:size(D,1)
            if (abs(dd2[i]-D[j])<R)
                dd2[i]=Inf; #Discard all Ritz values within a particular radius R
            end
        end
    end

    ii = sortperm(abs.(dd2-σ));

    mm_min = min(mm,length(ii));
    nu = dd2[ii[1:mm_min]];
    y = vv[:,ii[1:mm_min]];

    return nu,y;
end


# Residual-based Ritz value sorter:
# First discard all Ritz values within a distance R of any of the converged eigenvalues(of the original problem).
# Then select that Ritz value which gives the mm-th minimum product of (residual and distance from pole).
function residual_eigval_sorter(nep::NEP,dd,vv,σ,D,mm,R,Vk,errmeasure::Function=default_errmeasure(nep))

    eig_res = zeros(size(dd,1));
    dd2=copy(dd);

    ## Check distance of each eigenvalue of the projected NEP(i.e. in dd)
    ## from each eigenvalue that as already converged(i.e. in D)
    for i=1:size(dd,1)
        for j=1:size(D,1)
            if (abs(dd2[i]-D[j])<R)
                dd2[i]=Inf; #Discard all Ritz values within a particular radius R
            end
        end
    end

    #Compute residuals for each Ritz value
    for i=1:size(dd,1)
        eig_res[i] = errmeasure(dd[i],Vk*vv[:,i]);
    end

    #Sort according to methods
    ii = sortperm(eig_res .* abs.(dd2 .- σ))

    mm_min = min(mm,length(ii));

    nu = dd[ii[1:mm_min]];
    y = vv[:,ii[1:mm_min]];

    return nu,y;
end

# Threshold residual based Ritz value sorter:
# Same as residual_eigval_sorter() except that errors above a certain threshold are set to the threshold.
function threshold_eigval_sorter(nep::NEP,dd,vv,σ,D,mm,R,Vk,errmeasure::Function=default_errmeasure(nep),threshold=0.1)

    eig_res = zeros(size(dd,1));
    dd2=copy(dd);

    ## Check distance of each eigenvalue of the projected NEP(i.e. in dd)
    ## from each eigenvalue that as already converged(i.e. in D)
    for i=1:size(dd,1)
        for j=1:size(D,1)
            if (abs(dd2[i]-D[j])<R)
                dd2[i]=Inf; #Discard all Ritz values within a particular radius R
            end
        end
    end

    #Compute residuals for each Ritz value
    temp_res = 0;
    for i=1:size(dd,1)
        temp_res = errmeasure(dd[i],Vk*vv[:,i]);
        if(temp_res > threshold)
            eig_res[i] = threshold;
        else
            eig_res[i] = temp_res;
        end
    end

    #Sort according to methods
    ii = sortperm(eig_res.*abs.(dd2-σ));

    mm_min = min(mm,length(ii));

    nu = dd[ii[1:mm_min]];
    y = vv[:,ii[1:mm_min]];

    return nu,y;
end

nlar(nep::NEP;params...) = nlar(ComplexF64,nep::NEP;params...)
function nlar(::Type{T},
            nep::ProjectableNEP;
            orthmethod::Type{T_orth} = ModifiedGramSchmidt,
            nev::Int=10,                                     #Number of eigenvalues required
            submax_rest::Int=30,                             #Maximum dimension of subspace before restarting
            errmeasure::Function = default_errmeasure(nep),
            tol = eps(real(T))*100,
            maxit::Int = 100,
            λ0 = zero(T),
            v0 = randn(T,size(nep,1)),
            displaylevel::Int = 0,
            linsolvercreator::Function = default_linsolvercreator,
            R = 0.01,
            mm::Int = 4,
            eigval_sorter::Function = residual_eigval_sorter, #Function to sort eigenvalues of the projected NEP
            qrfact_orth::Bool = false,
            inner_solver_method = NEPSolver.DefaultInnerSolver) where {T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod}

        local σ::T = T(λ0); #Initial pole

        if (maxit > size(nep,1))
            @warn "Maximum iteration count maxit=$maxit larger than problem size n=$(size(nep,1)). Reducing maxit."
            maxit = size(nep,1);
        end

        λ0::T = T(λ0);
        #Initialize the basis V_1
        V = zeros(T, size(nep,1) ,maxit);
        X = zeros(T, size(nep,1) ,nev);
        V[:,1] = normalize(ones(T,size(nep,1)));
        Vk = zeros(T, size(nep,1) ,1);
        Vk[:,1] = V[:,1];

        D = zeros(T,nev);#To store the converged eigenvalues
        err_hyst=eps()*ones(maxit,nev) # Giampaolo's edit

        m = 0; #Number of converged eigenvalues
        k = 1;

        proj_nep = create_proj_NEP(nep);

        local linsolver::LinSolver = linsolvercreator(nep,σ);

        err = Inf;
        nu = λ0;
        u = v0;

        if(displaylevel == 1)
            println("##### Using inner solver:",inner_solver_method," #####");
        end

        while (m < nev) && (k < maxit)
            # Construct and solve the small projected PEP projected problem (V^H)T(λ)Vx = 0
            set_projectmatrices!(proj_nep,Vk,Vk);

            #Use inner_solve() to solve the smaller projected problem
            dd,vv = inner_solve(inner_solver_method,T,proj_nep,Neig=nev,σ=σ);

            # Sort the eigenvalues of the projected problem by measuring the distance from the eigenvalues,
            # in D and exclude all eigenvalues that lie within a unit disk of radius R from one of the
            # eigenvalues in D.
            nuv,yv = eigval_sorter(nep,dd,vv,σ,D,mm,R,Vk)
            nu = nuv[1]; y=yv[:,1];

            if (isinf(nu))
                error("We did not find any (non-converged) eigenvalues to target")
            end


            #Determine ritz vector and residual
            u = Vk*y; # Note: y and u are vectors (not matrices)

            #Normalize and compute residual
            normalize!(u);
            res = compute_Mlincomb(nep,nu,u);


            #Check for convergence of one of the eigenvalues
            err = errmeasure(nu,u);
            println(k," Error:",err," Eigval :",nu)
            err_hyst[k,m+1]=err;    # Giampaolo's edit
            if(err < tol)
                if(displaylevel == 1)
                    println("\n****** ",m+1,"th converged to eigenvalue: ",nu," errmeasure:",err,"  ******\n")
                end

                #Add to the set of converged eigenvalues and eigenvectors
                D[m+1] = nu;
                X[:,m+1] = u;

                ## Sort and select he eigenvalues of the projected problem as described before
                nuv,yv = eigval_sorter(nep,dd,vv,σ,D,mm,R,Vk)
                nu1=nuv[1];
                y1=yv[:,1];

                #Compute residual again
                u1 = Vk*y1;
                normalize!(u1);
                res = compute_Mlincomb(nep,nu1,u1);

                m = m+1;
            end

            #Compute new vector Δv to add to the search space V(k+1) = (Vk,Δv)
            Δv=lin_solve(linsolver,res)

            #Orthogonalize and normalize
            if (qrfact_orth)
                # Orthogonalize the entire basis matrix
                # together with Δv using QR-method.
                # Slow but robust.
                Q,_ = qr(hcat(Vk,Δv))
                Q = Matrix(Q)
                Vk=Q
                V[:,1:k+1]=Q;
                #println("Dist normalization:",opnorm(Vk'*Vk-I))
                #println("Size:",size(Vk), " N: ",opnorm(Vk[:,k+1]), " d:",opnorm(Δv))
            else
                h=zeros(T,k);
                orthogonalize_and_normalize!(Vk,Δv,h,orthmethod);

                #Expand basis
                V[:,k+1] = Δv;
                Vk = view(V,:,1:k+1);
            end

            #Check orthogonalization
            if(k < 100)
               println("CHECKING BASIS ORTHOGONALITY  ......     $(opnorm(Vk'*Vk - I))")
               #println("CHECKING ORTHO  ......     ",opnorm(Δv)," ....",h," .... ",g)
            end
            k = k+1;
        end


        #Throw no convergence exception if enough eigenvalues were not found.
        if(k >= maxit && m < nev)
            msg="Number of iterations exceeded. maxit=$(maxit) and only $(m) eigenvalues converged out of $(nev)."
            throw(NoConvergenceException(nu,u,err,msg))
        end

        return D,X,err_hyst;
    end
