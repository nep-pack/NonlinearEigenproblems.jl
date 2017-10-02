    export blocknewton

"""
    blocknewton

Applies the block Newton method to nep.

# References
* D. Kressner A block Newton method for nonlinear eigenvalue problems, Numer. Math., 114 (2) (2009), pp. 355-372
"""
function blocknewton(nep::AbstractSPMF;
                     S=zeros(3,3),
                     X=eye(size(nep,1),3),
#                       errmeasure::Function =                      default_errmeasure(nep::NEP),
                       tol=eps(real(eltype(S)))*100,
                       maxit=10,
                       displaylevel=0,
                       armijo_factor=1,
                       armijo_max=5)
    T=complex(eltype(S))
    # This implementation is for complex arithmetic only
    S=complex(S);
    X=complex(X);

    
    n=size(nep,1);
    p=size(S,1);

    RV_zeros=zeros(T,p,p);
    W=Vl(X,S);
    # reshape W to WW (3D matrix)
    WW=zeros(T,n,p,p)
    for j=1:size(W,2)
        WW[:,:,j]=W[(j-1)*n+(1:n),:];
    end

    # Main loop
    for k=1:maxit
        Res=  compute_MM(nep,S,X)
        resnorm=norm(Res);
        @ifd(@printf("Iteration %d: Resnorm: %e\n",k,resnorm))        
        if (resnorm<tol)
            return S,X
        end

        # Solve the linear system by first transforming S to
        # upper triangular form and then doing a backward substitution
        # and finally reversing the backward substitution 
        RR,QQ=schur(complex(S))
        dSt,dXt = newtonstep_linsys(nep, RR, X*QQ, WW, Res*QQ,
                                    0*QQ,
#                                    RV_zeros*QQ,
                                    displaylevel);
        dX=dXt*QQ';
        dS=QQ*dSt*QQ';


        # Update the approximations
        Xt=X-dX;
        St=S-dS;

        # Carry out the orthogonalization 
        (W,R)=qr(Vl(Xt,St),thin=true); 
        # reshape W to WW (3D matrix)
        for j=1:size(W,2)
            WW[:,:,j]=W[(j-1)*n+(1:n),:];
        end        
        X[:]=Xt/R; S[:]=(R*St)/R;
        
    end

    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(S,V,err,msg))

end



function newtonstep_linsys(nep::AbstractSPMF,S, X, W, RT, RV, displaylevel)
    T=complex(eltype(S))

    n = size(nep,1);
    p = size(X,2); l = size(W,3);
    dX = zeros(T, n,p); dS = zeros(T,p,p);

    # Fetch the SPMF 
    fv=get_fv(nep);
    Av=get_Av(nep)

    m = size(fv,1);

    # Initialize a tensor with a f_i(S) 
    fS = zeros(T,p,p,m); 
    for j = 1:m
        fS[:,:,j] = fv[j](S);
    end


    # Work column by-column
    for i = 1:p,
        s = S[i,i];
        T11=compute_Mder(nep,s)
        
        S_expanded=[S eye(S);zeros(S) s*eye(S)]
        # T12 = compute_MM(nep,[0*X X],S_expanded) # Can maybe be computed like this? Would avoid the explicit use of Av and and fv

        T12 = zeros(T,n,p);
        for j = 1:m
            DF = fv[j](S_expanded)
            DF1=DF[1:p,p+1:2*p];
            T12 = T12 + Av[j]*X*DF1;
        end
        T21 = W[:,:,1]'; 
        for j = 2:l
            T21 = T21 + s^(j-1) * W[:,:,j]'; 
        end
        DS = eye(p);
        T22 = zeros(T,p,p);
        for j = 2:l
            T22 = T22 + W[:,:,j]'*X*DS;
            DS = s*DS + S^(j-2) 
        end
        TT=[T11 T12; T21 T22]; # This can maybe be optimized with Schur complement
        sol =  TT \ [RT[:,i];RV[:,i]];
        # compute the dS and dX
        dX[:,i] = sol[1:n];
        dS[:,i] = sol[n+1:end];
        
        # Update RHS 
        Z = zeros(T,p,p);
        Z[:,i] = dS[:,i]; DS = Z;
        S2_expanded=[S Z;zeros(S) S];
        for j = 1:m
            
            DF = fv[j](S2_expanded)

            Z1a=dX[:,i]*(fS[i,i+1:p,j]).'
            Z1b=X*DF[1:p,p+i+1:2*p];
            Z1=(Z1a + Z1b );
            RT[:,i+1:p] =  RT[:,i+1:p]  - Av[j] * Z1;
        end
        for j = 2:l
            Z1a=dX[:,i]*  ((S^(j-2))[i,i+1:p].');
            Z1b=X*DS[:,i+1:p] ;
            Z1=W[:,:,j]' * ( Z1a +Z1b );
            RV[:,i+1:p] = RV[:,i+1:p] - Z1
            DS = DS*S + S^(j-2)*DS;
        end
    end
    return dS,dX
end


# Construct the orthogonalization matrix
# V=[X;X*S;XS^2;...]
function Vl(X,S)
    p=size(S,1); n=size(X,1);
    V=zeros(eltype(S),n*p,p);
    for j=1:p
        V[(j-1)*n+(1:n),:]=X*(S^(j-1));
    end
    return V
end

