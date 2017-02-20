#The non-linear Arnoldi method, as introduced in "An Arnoldi method for non-linear eigenvalue problems" by H.Voss

export nlar

###########################################################################################################

    function nlar(nep::NEP;
                nev=1,#Number of eigenvalues required
                errmeasure::Function =
                default_errmeasure(nep::NEP),
                tol=eps()*100,
                maxit=30,
                λ=0,
                v=randn(nep.n),
                displaylevel=0,
                nl_eigsolvertype=Union{AbstractString,Function},
                linsolvertype::DataType=DefaultLinSolver)

        σ = λ; #Initial pole 

        #Initialize the basis V_1
        V = zeros(nep.n,maxit);
        X = zeros(nep.n,nev);
        V[:,1] = ones(nep.n,1);
        V[:,1] = normalize(V[:,1]);
        Vk = V[:,1];

        D = zeros(nev,1);#To store the converged eigenvalues

        m = 0;#Number of converged eigenvalues
        k = 1;
        iter = 0;
        num_t = size(nep.A)[1]; #Number of monomial coefficients in the PEP = degree(PEP)+1
        while m < nev && iter < maxit 
            ### Construct the small projected PEP projected problem (V^H)T(λ)Vx = 0 using nl_eigsolvertype....(Currently works
            ### only for PEP) #####
            
            AA = Array{typeof(zeros(k,k))}(num_t);            
            for i=1:num_t
                AA[i] = zeros(k,k);
                if(k != 1)
                    AA[i] = Vk'*nep.A[i]*Vk;
                else
                    AA[i][:] = Vk'*nep.A[i]*Vk;
                end
            end

            pep_proj = PEP(AA);

            #Solve the projected problem
            ν,y =res_inv(pep_proj,maxit=30,displaylevel=1);
            print("\n\n")
        
            if(k == 1)
                y = y[1];
            end
            #Determine ritz vector and residual
            u = Vk*y; 
            res = compute_Mder(nep,ν,0)*u;

            #Check for convergence of one of the eigenvalues
            err = errmeasure(ν,u);
            if(err < tol)
                if(displaylevel == 1)
                    println("Eigenvalue: ",λ," errmeasure:",err)
                end
                D[m+1] = ν;
                X[:,m+1] = u;
                m = m+1; 
                k = 0;
            end

            #Compute new vector Δv to add to the search space V(k+1) = (Vk,Δv)
            local linsolver::LinSolver=linsolvertype(compute_Mder(nep,σ));
            Δv=lin_solve(linsolver,res,tol=1e-16)
            V[:,k+1] = Δv;

            #Orthogonalize
            if(k != 0)
                h = V[:,1:k]'*Δv;
                V_k1 = Δv-V[:,1:k]*h;
                V[:,k+1] = V_k1;
            end

            #Expand
            Vk = V[:,1:k+1];

            k = k+1; 

            iter = iter+1;
        end

        return D,X;
    end

