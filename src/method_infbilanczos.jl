
export infbilanczos
"""
    The Infinite Bi-Lanczos
"""
    infbilanczos(nep::NEP,nept::NEP;params...)=infbilanczos(Complex128,nep,nept;params...)
    function infbilanczos{T<:Number}(::Type{T},
                          nep::NEP,
                          nept::NEP;  # Transposed NEP
                          maxit=30,
                          linsolvercreator::Function=default_linsolvercreator,
                          linsolvertcreator::Function=linsolvercreator,
                          v=randn(real(T),size(nep,1)),
                          u=randn(real(T),size(nep,1)),
                          tol=1e-12,
                          Neig=maxit,                                  
                          errmeasure::Function = default_errmeasure(nep::NEP),
                          σ=0.0,
                          γ=1,
                          displaylevel=0)

        
        n=size(nep,1);
        σ=T(σ);

        # Linear systems solver for both M(σ) and M(σ)^H

        # Shift σ \neq 0 not implemented
        
        local M0inv::LinSolver = linsolvercreator(nep,σ);
        local M0Tinv::LinSolver = linsolvertcreator(nept,conj(σ));        
        
        #        
        m=maxit;
        qt=lin_solve(M0Tinv,u);
        q=v;
        q=q/dot(qt,compute_Mlincomb(nep,σ,q,ones(1),1));
        
        Q0=zeros(T,n,m);                  # represents Q_{k-1}
        Qt0=zeros(T,n,m);                 # represents \til{Q}_{k-1}
        R1=zeros(T,n,m); R1[:,1]=q;       # represents R_{k}
        Rt1=zeros(T,n,m); Rt1[:,1]=qt;    # represents \tild{R}_{k}  
        Z2=zeros(T,n,m);
        Zt2=zeros(T,n,m); 
        Q_basis=zeros(T,n,m);
        Qt_basis=zeros(T,n,m);


        R2=zeros(T,n,m); # Needed?
        Rt2=zeros(T,n,m); # Needed?        
        
        Q1=zeros(T,n,m); # Needed?
        Qt1=zeros(T,n,m); # Needed?


        # Vectors storing the diagonals
        alpha=zeros(T,m);               
        beta=zeros(T,m);
        gamma=zeros(T,m);


        @ifd(@printf("%e %e\n",norm(q), norm(qt)));

        k=1;

        @ifd(@printf("Iteration:"));
        while k < m
            @ifd(@printf("%d ", k));
            # Note: conjugate required since we compute s'*r not r'*s
            omega = conj(left_right_scalar_prod(T,nep,nept,Rt1,R1,k,k,σ));

            beta[k] = sqrt(abs(omega));        
            gamma[k] = conj(omega) / beta[k];
            
            # Step 11-12
            
            Q1[:,1:k]=R1[:,1:k]/beta[k];
            Qt1[:,1:k]=Rt1[:,1:k]/conj(gamma[k]);

            #@printf("\n%e %e\n", Q1[1,1],Qt1[1,1]);
            
            
            # Extra step, to compute Ritz vectors eventually
            Q_basis[:,k] = Q1[:,1];
            Qt_basis[:,k] = Qt1[:,1];        

             # Step 1: Compute Z_{k+1} 
            Dk=diagm(1./(exp.(lfact.(1:k))));
            b1_tmp=compute_Mlincomb(nep,σ,Q1[:,1:k]*Dk,ones(k),1);
            b1=-lin_solve(M0inv,b1_tmp);
            Z2[:,k] = b1;
            


            # Step 2: Compute \til{Z}_{k+1} 
            bt1_tmp=compute_Mlincomb(nept,conj(σ),Qt1[:,1:k]*Dk,ones(k),1);
            bt1=-lin_solve(M0Tinv,bt1_tmp);
            Zt2[:,k] = bt1

            # Step 3: Compute R_{k+1}            
            R2[:,1] = Z2[:,k];            
            R2[:,2:(k+1)]=Q1[:,1:k];
            if k > 1
                R2[:,1:(k-1)]=R2[:,1:(k-1)]-gamma[k]*Q0[:,1:(k-1)];
            end

            # Step 4: Compute \til{R}_{k+1}
            Rt2[:,1] = Zt2[:,k];
            Rt2[:,2:(k+1)]=Qt1[:,1:k];
            if k > 1
                Rt2[:,1:(k-1)]=Rt2[:,1:(k-1)]-conj(beta[k])*Qt0[:,1:(k-1)];
            end

            # Step 5: Compute \alpha_k
            alpha[k+1]=left_right_scalar_prod(T,nep,nept,Qt1,R2,k,k+1,σ);
            

            #Step 6: Compute R_{k+1}
            R2[:,1:k]=R2[:,1:k]-alpha[k+1]*Q1[:,1:k];


 
            #Step 7: Compute \til{R}_{k+1}
            Rt2[:,1:k]=Rt2[:,1:k]-conj(alpha[k+1])*Qt1[:,1:k];
            

            # shift  the matrices: 
            (R1,R2)=(R2,R1);  R2[:,:]=0  # swap and reset forgotten variable
            (Rt1,Rt2)=(Rt2,Rt1);  Rt2[:,:]=0
            (Q0,Q1)=(Q1,Q0);  Q1[:,:]=0
            (Qt0,Qt1)=(Qt1,Qt0);  Qt1[:,:]=0

            
            k=k+1;
           
        end
        omega = left_right_scalar_prod(T,nep,nept,Rt1,R1,m,m,σ);
        
        beta[m] = sqrt(abs(omega));        
        gamma[m] = conj(omega) / beta[m];
        alpha=alpha[2:end];  # \alpha_1 stored in alpha(2)
        beta=beta[2:end];    # we do not need \beta_1
        gamma=gamma[2:end];  # we do not need \gamma_1
        

        @ifd(@printf("done \n"));
        
        TT = full(spdiagm((beta[1:m-1],alpha[1:m-1],gamma[1:m-1]), -1:1));
        λ = σ+1./eigvals(TT);
        return λ,0,TT;
    end

    function left_right_scalar_prod{T}(::Type{T}, nep,nept,At,B,ma,mb,σ)
        # Compute the scalar product based on the function nep.M_lin_comb  
        c=0;
        # This is the nasty double loop, which brings
        # complexity O(m^3n). Will be limiting if we do many iterations
        XX=zeros(T,size(B,1),mb); # pre-allocate
        for j=1:ma
            #dd=1./factorial(j:(j+mb-1));
            dd=1./exp.(lfact.(j:(j+mb-1)));
            XX=broadcast(*,B[:,1:mb],dd'); # diag scaling
            #XX=bsxfun(@times,B(:,1:mb),dd);  # Column scaling: Faster than constructing

            # compute Mlincomb starting from derivative j
            z=-compute_Mlincomb(nep,σ,XX,ones(size(XX,2)),j);
            c=c+dot(At[:,j],z);
        end  
        return c
    end

