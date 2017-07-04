#The non-linear Arnoldi method, as introduced in "An Arnoldi method for non-linear eigenvalue problems" by H.Voss

export nlar

################################################################################################################

function default_proj_solver(pnep::Proj_NEP,nev=10,σ=0.0)
    if (isa(pnep,Proj_PEP))   
        return polyeig(pnep.nep_proj)
    else
        λ,Q,err=iar(pnep,Neig=2*nev+3,σ=σ,maxit=100)
        return λ,Q
    end
end

## D = already computed eigenvalues
## dd, vv eigenpairs of projected problem
## σ target
## returns at most mm values
function  default_eigval_sorter(dd,vv,σ,D,mm)
    R=0.01;
    dd2=copy(dd);
    for i=1:size(dd,1)
        for j=1:size(D,1)
            if (abs(dd2[i]-D[j])<R)
                dd2[i]=Inf;
            end
        end
    end            
    ii = sortperm(abs(dd2-σ));

    mm_min=min(mm,length(ii));
    nu = dd2[ii[1:mm_min]];
    y = vv[:,ii[1:mm_min]];
    
    return nu,y
end

function nlar(nep::NEP;
                nev=10,#Number of eigenvalues required
                errmeasure::Function =
                default_errmeasure(nep),
                tol=1e-6,
                maxit=100,
                λ=0,
                v=randn(size(nep,1)),
                displaylevel=0,
                nl_eigsolvertype=Union{AbstractString,Function},
                  linsolvercreator::Function=default_linsolvercreator,
                  proj_solve::Function = default_proj_solver,
              eigval_sorter::Function = default_eigval_sorter,
              qrfact_orth::Bool=false)

        local σ::Complex128 = λ; #Initial pole 

        if (maxit>size(nep,1))
            warn("Maximum iteration count maxit="*string(maxit)*" larger than problem size n="*string(size(nep,1))*". Reducing maxit.")
            maxit=size(nep,1);
        end
        #Initialize the basis V_1
        V = zeros(Complex128, size(nep,1) ,maxit);
        X = zeros(Complex128, size(nep,1) ,nev);
        V[:,1] = normalize(ones(size(nep,1)));
        Vk = zeros(Complex128, size(nep,1) ,1);
        Vk[:,1] = V[:,1];

        D = zeros(Complex128,nev);#To store the converged eigenvalues

        m = 0; #Number of converged eigenvalues
        k = 1; 

        proj_nep=create_proj_NEP(nep);
        
        local linsolver::LinSolver=linsolvercreator(nep,σ);
 

        ### TODO: What happens when k reaches maxit? NoConvergenceError? ###
        while (m < nev) && (k < maxit)
            ### Construct the small projected PEP projected problem (V^H)T(λ)Vx = 0 using nl_eigsolvertype....(Currently works
            ### only for PEP) #####

            set_projectmatrices!(proj_nep,Vk,Vk);
            
            dd,vv = proj_solve(proj_nep,nev,σ);            

            nuv,yv = eigval_sorter(dd,vv,σ,D, 4)
            

            nu=nuv[1];
            y=yv[:,1];

            if (isinf(nu))
                error("We did not find any (non-converged) eigenvalues to target")
            end
            

            #Determine ritz vector and residual
            u = Vk*y; # Note: y and u are vectors (not matrices)

            u = normalize(u);
            res = compute_Mlincomb(nep,nu,u);

           
            #Check for convergence of one of the eigenvalues
            err = errmeasure(nu,u);
            println(k," Error:",err," Eigval :",nu)
            if(err < tol)
                if(displaylevel == 1)
                    println("\n\n****** ",m+1,"th converged to eigenvalue: ",nu," errmeasure:",err,"  ******\n")
                end
                D[m+1] = nu;
                X[:,m+1] = u;

                #Change the pole
                #σ = 2.0*ν;

                nuv,yv = eigval_sorter(dd,vv,σ,D, 4)
                nu1=nuv[1];
                y1=yv[:,1];
                
                #Compute residual again
                u1 = Vk*y1; 
                u1 = normalize(u1);
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
                Q,R=qr(hcat(Vk,Δv),thin=true)
                Vk=Q
                V[:,1:k+1]=Q;
                #println("Dist normalization:",norm(Vk'*Vk-eye(k+1)))
                #println("Size:",size(Vk), " N: ",norm(Vk[:,k+1]), " d:",norm(Δv))
            else

                # Do our own (double) Gram-Schmidt
                h = V[:,1:k]'*Δv;
                V_k1 = Δv-V[:,1:k]*h;
                g = V[:,1:k]'*V_k1;
                V_k1 = V_k1-V[:,1:k]*g;
                V_k1 = normalize(V_k1)
               
                #Expand
                V[:,k+1] = V_k1;
                Vk = V[:,1:k+1];
            end

            #Check orthogonalization
            if(k < 100)
               #println("CHECKING ORTHO  ......     ",norm(Vk'*Vk-eye(Complex128,k+1)),"\n\n")
                #println("CHECKING ORTHO  ......     ",norm(Δv)," ....",h," .... ",g,"\n") 
            end
            k = k+1; 
        end

        return D,X;
    end

