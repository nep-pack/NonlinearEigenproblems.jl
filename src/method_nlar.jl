#The non-linear Arnoldi method, as introduced in "An Arnoldi method for non-linear eigenvalue problems" by H.Voss

export nlar

################################################################################################################

    function nlar(nep::NEP;
                nev=10,#Number of eigenvalues required
                errmeasure::Function =
                default_errmeasure(nep::NEP),
                tol=1e-6,
                maxit=100,
                λ=0,
                v=randn(nep.n),
                displaylevel=0,
                nl_eigsolvertype=Union{AbstractString,Function},
                linsolvertype::DataType=DefaultLinSolver)

        σ = λ; #Initial pole 

        if (maxit>size(nep,1))
            warn("Maximum iteration count maxit="*string(maxit)*" larger than problem size n="*string(size(nep,1))*". Reducing maxit.")
            maxit=size(nep,1);
        end
        #Initialize the basis V_1
        V = zeros(Complex128,nep.n,maxit);
        X = zeros(Complex128,nep.n,nev);
        V[:,1] = normalize(ones(nep.n));
        Vk = V[:,1];

        D = zeros(Complex128,nev);#To store the converged eigenvalues

        m = 0;#Number of converged eigenvalues
        k = 1;


        qrmethod_orth=true;  # 
        
        local linsolver::LinSolver=linsolvertype(compute_Mder(nep,σ,0));
 
        num_t = size(nep.A)[1]; #Number of monomial coefficients in the PEP = degree(PEP)+1
        while m < nev 
            ### Construct the small projected PEP projected problem (V^H)T(λ)Vx = 0 using nl_eigsolvertype....(Currently works
            ### only for PEP) #####
            
            AA = Array{typeof(zeros(Complex128,k,k))}(num_t);            
            for i=1:num_t
                AA[i] = zeros(Complex128,k,k);
                if(k != 1)
                    AA[i] = Vk'*nep.A[i]*Vk;
                else
                    AA[i][:] = Vk'*nep.A[i]*Vk;
                end
            end

            pep_proj = PEP(AA);

            #Solve the projected problem
            dd,vv = polyeig(pep_proj);

            ii = sortperm(abs(dd-σ));

            ν = dd[ii[m+1]];
            y = vv[:,ii[m+1]];

            if(k == 1)
                y = y[1]
            end
            
            #Determine ritz vector and residual
            u = Vk*y; 

            u = normalize(u);
            res = compute_Mlincomb(nep,ν,u);

           
            #Check for convergence of one of the eigenvalues
            err = errmeasure(ν,u);
            println(k," Error:",err," Eigval :",ν)
            if(err < tol)
                if(displaylevel == 1)
                    println("\n\n****** ",m+1,"th converged to eigenvalue: ",ν," errmeasure:",err,"  ******\n")
                end
                D[m+1] = ν;
                X[:,m+1] = u;

                #Change the pole
                #σ = 2.0*ν;

                #Compute residual again
                ν1 = dd[ii[m+2]];
                y1 = vv[:,ii[m+2]]
                u1 = Vk*y1; 
                u1 = normalize(u1);
                res = compute_Mlincomb(nep,ν1,u1);

                m = m+1; 
            end

            #Compute new vector Δv to add to the search space V(k+1) = (Vk,Δv)

            Δv=lin_solve(linsolver,res)

            #Orthogonalize and normalize
            if (qrmethod_orth)

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

