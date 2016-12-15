function result=iar(nep,startvec,MAXN)
%   A reference implementation of the infinite Arnoldi method
%   taylor version. 
% 
    n=length(startvec);

    
    mu=0; % In Bi-Lanczos comparison we only use mu=0 since only
          % implemented in nep
    
    
    x0=startvec; x0=x0/norm(x0);
    H=zeros(MAXN+1,MAXN);
    fprintf('IAR Iteration:');
    V=zeros(n*(MAXN+1),MAXN);
    V(1:n,1)=x0; 
    t0=tic();
    result.timing_iteration=[];
    for k=1:MAXN 
        fprintf('%d ', k);

        % Create L matrix: Equation (17) in nummath paper
        L=diag(1./(1:k));
        vecx=V(1:(k*n),k);
        X=reshape(vecx,n,k);
        
        % Compute y1,..yn: Equation (13) in nummath paper
        Y=zeros(size(X,1),size(X,2)+1);
        Y(:,2:end)=X*L;
        
        % Compute y0: Equation (18) in nummath paper
        % y0=nep.M_lin_comb(X,Y,mu);
        z0=nep.M_lin_comb(Y(:,2:end),1);        
        y0=-nep.M_solve(z0);
        Y(:,1)=y0;
        vecY=reshape(Y,n*(k+1),1);

        % Expand V to prepare for orthogonalization
        %V=[V;zeros(n,N)];   % Needs to be pre-allocated in CPU-comparison
        
        
        % Orthogonalize 
        V0=V(1:((k+1)*n),1:k);
        h=V0'*vecY; 
        vecY=vecY-V0*h;
 
        g=V0'*vecY;   % twice to be sure
        vecY=vecY-V0*g;        
        h=h+g;
        
        beta=norm(vecY);  % normalize
        vecY=vecY/beta;
        H(1:(k+1),k)=[h;beta];
        V(1:((k+1)*n),k+1)=vecY;
        result.timing_iteration=[result.timing_iteration;toc(t0)];
    end
    fprintf('  done. \n');
    
    
    result.V=V;
    result.H=H;

    [VV,DD]=eig(result.H(1:end-1,1:end));;
    ritzvals= diag(DD);
    
    [Y,I]=sortrows(-abs(ritzvals)); ritzvals=ritzvals(I);
    
    result.eigvals=mu+1./ritzvals;

    X=V(1:n,1:MAXN)*VV;
    X=X(:,I);
    for i=1:size(X,2)
        X(:,i)=X(:,i)/norm(X(:,i));
    end
    result.X=X;

    if (0)  % for debugging only
       result.Xfull=V(:,1:MAXN)*VV; 
    end
   