using Combinatorics

export infbilanczos
"""
    The Infinite Bi-Lanczos
"""
    function infbilanczos(nep::NEP,
                          nept::NEP;  # Transposed NEP
                          maxit=30,
                          linsolvercreator::Function=default_linsolvercreator,
                          tol=1e-12,
                          Neig=maxit,                                  
                          errmeasure::Function = default_errmeasure(nep::NEP),
                          σ=0.0,
                          γ=1,
                          displaylevel=0)

        
        n=size(nep,1);
#        
        m=maxit;
        q=ones(n);
        qt=ones(n);
        Q0=zeros(n,m);                  # represents Q_{k-1}
        Qt0=zeros(n,m);                 # represents \til{Q}_{k-1}
        R1=zeros(n,m); R1[:,1]=q;       # represents R_{k}
        Rt1=zeros(n,m); Rt1[:,1]=qt;    # represents \tild{R}_{k}  
        Z2=zeros(n,m);
        Zt2=zeros(n,m); 
        Q_basis=zeros(n,m);
        Qt_basis=zeros(n,m);

        alpha=zeros(m);               
        beta=zeros(m);
        gamma=zeros(m);
        itertime=zeros(m);



        

        return 0,0,0;
    end

    function left_right_scalar_prod(nep,nept,At,B,ma,mb)
        # Compute the scalar product based on the function nep.M_lin_comb  
        c=0;
        # This is the nasty double loop, which brings
        # complexity O(m^3n). Will be limiting if we do many iterations
        XX=zeros(size(B,1),mb); # pre-allocate
        for j=1:ma
            #dd=1./factorial(j:(j+mb-1));
            dd=1./exp(lfact(j:(j+mb-1)));
            XX=broadcast(*,B[:,1:mb],dd'); # diag scaling
            #XX=bsxfun(@times,B(:,1:mb),dd);  # Column scaling: Faster than constructing
            # compute Mlincomb starting from derivative j
            z=-compute_Mlincomb(nep,XX,ones(size(XX,2)),j);
            c=c+At[:,j]'*z;
        end  
        return c
    end

