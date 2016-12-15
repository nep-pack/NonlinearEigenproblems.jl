function result=tiar(compute_y0,mu,startvec,MAXN)
% The tensor infinite Arnoldi method (TIAR) as described in
% Algorithm 2 in http://arxiv.org/abs/1502.01613
%
% Syntax same as iar.m 
% The user must specify target mu and a function compute_y0(X,Y,mu)% 

    global timing_y0;
    
    n=length(startvec);

    x0=startvec;
    x0=x0/norm(x0);
    Z=zeros(n,(MAXN+1));
    Z(1:n,1)=x0;
    H=zeros(MAXN+1,MAXN);

    a=NaN*zeros(MAXN+1,MAXN+1,MAXN+1);
    a(1,1,1)=1;

    Y=NaN*zeros(n,MAXN+1); % Pre-allocate
    X=NaN*zeros(n,MAXN);  % Pre-allocate
    use_selective_reorth=0; % Modified 2016-03-01 
    if (use_selective_reorth)
         REORTH_TOL=eps*10;   % Warning: There seems to be
                              % round-off problems with selective reorth
    else
        REORTH_TOL=0; 
    end
    fprintf('TIAR Iteration:');
    timing.t0=0;
    timing.t1=0;   
    timing.t1_0=0;   
    timing.t1_0a=0;   
    timing.t1_0b=0;   
    timing.t1_0c=0;   
    timing.t1_1=0;   
    timing.t2=0;
    timing.t3=0;
    timing.t4=0;
    for k=1:MAXN 
        fprintf('%d ', k);

        tt0=tic;
        %% Step 3: Compute y_2... y_{k+1}        
        ZZ=Z(:,1:k);
        for j=2:k+1             
            z=zeros(n,1);
            aa=sparse(reshape(a(j-1,k,1:k),k,1)); 
            % use sparse to avoid 0 mults. Described as a
            % matrix-matrix multiply in paper. Not sure if this
            % is faster than actual matrix matrix multiplies in
            % practice. Probably depends on architecture.  
            
            z=ZZ*aa;
            
            Y(:,j)=(1/(j-1))*z;
            
            %X(:,j-1)=z; % Not necessary
        end
        timing.t0=timing.t0+toc(tt0);
            
        
        %% Step 3:  Compute y1  
        tt1=tic;
        y=compute_y0(X(:,1:k),Y(:,1:k+1),mu);
        timing.t1=timing.t1+toc(tt1);
        timing.t1_0=timing.t1_0+timing_y0.t0;
        timing.t1_1=timing.t1_1+timing_y0.t1;        
        timing.t1_0a=timing.t1_0a+timing_y0.t0a;
        timing.t1_0b=timing.t1_0b+timing_y0.t0b;        
        timing.t1_0c=timing.t1_0c+timing_y0.t0c;        
        
        %% Step 4: Orthogonalize y against Z (compute t and zk+1
        tt2=tic;
        yorg=y;
        
        t=Z(1:n,1:k)'*y;
        y=y-Z(1:n,1:k)*t;
        
        tt=Z(1:n,1:k)'*y;

        if (norm(tt)>REORTH_TOL);
            y=y-Z(1:n,1:k)*tt;
            t=t+tt;
        end
        b=norm(y);
        t=[t;b];
        Z(:,k+1)=y/b;
        timing.t2=timing.t2+toc(tt2);
        
        %% Step 5: Compute G 
        tt3=tic;
        G=NaN*zeros(k+1,k+1);
        % (4.8a)
        for l=1:k+1
           G(1,l)=t(l);    
        end
        % (4.8b)
        for i=2:(k+1)
            for l=1:k
                G(i,l)=(1/(i-1))*a(i-1,k,l);
            end
        end
        % (4.8c)
        for i=2:k+1
            G(i,k+1)=0;    
        end
        
        %% Step 6: Compute h
        
        h=zeros(k,1);
        for l=1:k
            h=h + a(1:k,1:k,l)'*G(1:k,l);   % (4.10)
        end
        timing.t3=timing.t3+toc(tt3);

        
        %% Step 7: Compute F
        tt4=tic;
        F=NaN*zeros(k+1,k+1);
        for l=1:k
            F(:,l)=G(:,l)-[a(1:k,1:k,l);zeros(1,k)]*h;
        end
        F(:,k+1)=G(:,k+1);
        
        
        %% Step 8: Possibly reorthogonalize 
        hh=zeros(k,1);
        for l=1:k
            hh=hh + a(1:k,1:k,l)'*F(1:k,l);   % (4.10)
        end
        if (norm(hh)>REORTH_TOL)
            FF=NaN*zeros(k+1,k+1);
            for l=1:k
                FF(:,l)=F(:,l)-[a(1:k,1:k,l);zeros(1,k)]*hh;
            end
            FF(:,k+1)=F(:,k+1);
            F=FF;
            h=h+hh;
        end
        
        %% Step 9: Compute beta
        
        beta=norm(F,'fro');

        %% Step 10: expand the tensor a
        
        for i=1:k+1
            for l=1:k+1
                a(i,k+1,l)=F(i,l)/beta;   % (4.14a)
                if (i<=k)
                    a(k+1,i,l)=0;         % (4.14b)
                    a(l,i,k+1)=0;         % (4.14c)
                end
            end
        end

        %% Step 12:
        H(1:(k+1),k)=[h;beta];
        timing.t4=timing.t4+toc(tt4);

    end
    fprintf('  done. \n');    
    fprintf(' timing:t0=%f  t1=%f t2=%f t3=%f t4=%f\n',timing.t0,timing.t1,timing.t2,timing.t3,timing.t4);
    fprintf(' timing:t1_0=%f t1_1=%f\n',timing.t1_0,timing.t1_1);
    fprintf(' timing:t1_0a=%f t1_0b=%f t1_0c=%f\n',...
            timing.t1_0a,timing.t1_0b,timing.t1_0c);
    result.H=H;
    
    
    %% Post-calculation to get the eigenvalue approximations
    
    [U,D]=eig(H(1:end-1,1:end));
    evps=mu+1./diag(D);  % reverse shift
     
    % Set eigval approximations
    result.eigvals=evps; 
    V1=zeros(n,k);

    for i=1:k
        V1(:,i)=Z*reshape(a(1,i,:),k+1,1);
    end
    result.V1=V1;

    % set eigvec approximations
    result.X=V1(:,1:k)*U;

    % normalize eigvec approximations
    for j=1:size(result.X,2)
       result.X(:,j)=result.X(:,j)/norm(result.X(:,j));
    end
