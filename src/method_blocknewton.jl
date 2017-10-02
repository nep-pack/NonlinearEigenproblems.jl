    export blocknewton


function blocknewton(nep::NEP;S=eye(3,3),X=eye(size(nep,1),3))

    T=Complex128
    n=size(nep,1);
    p=size(S,1);

    RV_zeros=zeros(T,p,p);
    W=Vl(X,S);
    # reshape W to WW (3D matrix)
    WW=zeros(T,n,p,p)
    for j=1:size(W,2)
        WW[:,:,j]=W[(j-1)*n+(1:n),:];
    end
    
    
    hist=[];

    resnorm=1;
    Shist=[];

    k=0;
    kmax=10;
    while ((resnorm > sqrt(eps())) && k<kmax)
        Res=  compute_MM(nep,S,X)
        
        RR,QQ=schur(complex(S))

        dSt,dXt = newtonstep_linsys(nep, RR, X*QQ, WW, Res*QQ, RV_zeros*QQ);
        dX=dXt*QQ';
        dS=QQ*dSt*QQ';
        Xt=X-dX;
        St=S-dS;

        (W,R)=qr(Vl(Xt,St),thin=true); 
        

        # reshape W to WW (3D matrix)
        for j=1:size(W,2)
            WW[:,:,j]=W[(j-1)*n+(1:n),:];
        end
        
        X=Xt*inv(R); S=(R*St)*inv(R);
        
        resnorm=norm(compute_MM(nep,S,X));
        #hist=[hist;resnorm];
        #Shist=[Shist;S];
        k=k+1;
        println("norm(S)=",norm(S))
        println("norm(X)=",norm(X))        
        println("*** RESNORM:",resnorm)
        
    end
    return S,X


end



function newtonstep_linsys(nep,S, X, W, RT, RV )
    T=Complex128

    n = size(nep,1);
    k = size(X,2); l = size(W,3);
    dX = zeros(T, n,k); dS = zeros(T,k,k);

    fv=get_fv(nep);
    m = size(fv,1);
    Av=get_Av(nep)
    fS = zeros(Complex128,k,k,m); 
    for j = 1:m
        fS[:,:,j] = fv[j](S);
    end


    # Work column by-column
    for i = 1:k,
        s = S[i,i];
        T11=compute_Mder(nep,s)
        
        S_expanded=[S eye(S);zeros(S) s*eye(S)]
        # T12 = compute_MM(nep,[0*X X],S_expanded) # Can maybe be computed like this?
        println("m=",m)

        T12 = zeros(T,n,k);
        for j = 1:m
            DF = fv[j](S_expanded)
            DF1=DF[1:k,k+1:2*k];
            T12 = T12 + Av[j]*X*DF1;
            
            println("T12:",T12, " DF1:",DF1, " X:",X)
        end
        T21 = W[:,:,1]'; 
        for j = 2:l
            T21 = T21 + s^(j-1) * W[:,:,j]'; 
        end
        DS = eye(k);
        T22 = zeros(T,k,k);
        for j = 2:l
            println("T22:",T22)
            T22 = T22 + W[:,:,j]'*X*DS;
            DS = s*DS + S^(j-2) # pS(:,:,j-1);
        end
        println("k:",k)        
        println("size T11:",size(T11))
        println("size T12:",size(T12))
        println("size T21:",size(T21))
        println("size T22:",size(T22))
        println("T22:",T22)
        TT=[T11 T12; T21 T22];
        println("TT:",TT)
        sol =  TT \ [RT[:,i];RV[:,i]];
        dX[:,i] = sol[1:n];
        dS[:,i] = sol[n+1:end];
        # Update right-hand side
        Z = zeros(T,k,k);
        Z[:,i] = dS[:,i]; DS = Z;
        S2_expanded=[S Z;zeros(S) S];
        for j = 1:m
            println("size(S2_expanded):",size(S2_expanded))
            
            DF = fv[j](S2_expanded)#feval(f,j,[S, Z;zeros(k) S ]);

            Z1a=dX[:,i]*(fS[i,i+1:k,j]).'
            println("typeof(DF):",typeof(DF))
            Z1b=X*DF[1:k,k+i+1:2*k];
            Z1=(Z1a + Z1b );
            println("size(Z1)=",size(Z1))
            RT[:,i+1:k] =  RT[:,i+1:k]  - Av[j] * Z1;
        end
        for j = 2:l

            #println("SS:",(((S^(j-2))[i,i+1:k]).'))
            Z1a=dX[:,i]*  ((S^(j-2))[i,i+1:k].');
            Z1b=X*DS[:,i+1:k] ;
            Z1=W[:,:,j]' * ( Z1a +Z1b );
            RV[:,i+1:k] = RV[:,i+1:k] - Z1
            DS = DS*S + S^(j-2)*DS;
        end
    end
    return dS,dX
end



function Vl(X,S)
    p=size(S,1); n=size(X,1);
    V=zeros(Complex128,n*p,p);
    for j=1:p
        V[(j-1)*n+(1:n),:]=X*(S^(j-1));
    end
    return V
end

