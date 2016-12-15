function result=infbilanczos(nep,q,qt,m)
%  Infinite Bi-Lanczos for nonlinear eigenvalue problems
%
%  This is an implementation of the infinite Bi-Lanczos procedure
%  described in 
% 
%  The infinite Bi-Lanczos method for nonlinear eigenvalue problems
%  Sarah W. Gaaf and Elias Jarlebring 
%  Arxiv prepring, 2016
% 
%
    
    k=1;
    n=length(q);
    Q0=zeros(n,m);                  % represents Q_{k-1}
    Qt0=zeros(n,m);                 % represents \til{Q}_{k-1}
    R1=zeros(n,m); R1(:,1)=q;       % represents R_{k}
    Rt1=zeros(n,m); Rt1(:,1)=qt;    % represents \tild{R}_{k}  
    Z2=zeros(n,m);
    Zt2=zeros(n,m); 
    Q_basis=zeros(n,m);
    Qt_basis=zeros(n,m);

    alpha=zeros(m,1);               
    beta=zeros(m,1);
    gamma=zeros(m,1);
    itertime=zeros(m,1);
    
    
    % Allow the possibility to supply a user specified 
    % left-right scalar product
    if (isfield(nep,'left_right_scalar_prod'))
        fprintf('Using user-specified left_right_scalar_prod\n');
        this_left_right_scalar_prod=nep.left_right_scalar_prod;
    else
        fprintf('Using generic left_right_scalar_prod\n');
        this_left_right_scalar_prod=@left_right_scalar_prod;
    end

    t0 = tic();
    result.timing_iteration=NaN*zeros(m-1,1);
    result.timing_scalarprod=NaN*zeros(m-1,1);
    fprintf('Infinite Bi-Lanczos iteration:');
    while k < m
        fprintf('%d ', k);
        % Step 8-10
        st0=tic;
        omega = conj(this_left_right_scalar_prod(nep,Rt1,R1,k,k)); % Note: conjugate required since we compute s'*r not r'*s
        result.timing_scalarprod(k)=toc(st0);
        beta(k) = sqrt(abs(omega));        
        gamma(k) = conj(omega) / beta(k);
        
        % Step 11-12
        Q1(:,1:k)=R1(:,1:k)/beta(k);
        Qt1(:,1:k)=Rt1(:,1:k)/conj(gamma(k));
        
        % Extra step, to compute Ritz vectors eventually
        Q_basis(:,k) = Q1(:,1);
        Qt_basis(:,k) = Qt1(:,1);        
        
        % Step 1: Compute Z_{k+1} 
        Dk=diag(1./(factorial(1:k)));
        b1_tmp=nep.M_lin_comb(Q1(:,1:k)*Dk,1);
        b1=-nep.M_solve(b1_tmp);
        Z2(:,k) = b1;
        
        % Step 2: Compute \til{Z}_{k+1} 
        bt1_tmp=nep.Mt_lin_comb(Qt1(:,1:k)*Dk,1);
        bt1=-nep.Mt_solve(bt1_tmp);
        Zt2(:,k) = bt1;
        
        % Step 3: Compute R_{k+1}
        R2(:,1) = Z2(:,k);
        R2(:,2:(k+1))=Q1(:,1:k);
        if k > 1; R2(:,1:(k-1))=R2(:,1:(k-1))-gamma(k)*Q0(:,1:(k-1)); end
        
        % Step 4: Compute \til{R}_{k+1}
        Rt2(:,1) = Zt2(:,k);
        Rt2(:,2:(k+1))=Qt1(:,1:k);
        if k > 1; Rt2(:,1:(k-1))=Rt2(:,1:(k-1))-conj(beta(k))*Qt0(:,1:(k-1)); end
        
        % Step 5: Compute \alpha_k
        st0=tic();
        alpha(k+1)=this_left_right_scalar_prod(nep,Qt1,R2,k,k+1);
        result.timing_scalarprod(k)=result.timing_scalarprod(k)+toc(st0);
        
        %Step 6: Compute R_{k+1}
        R2(:,1:k)=R2(:,1:k)-alpha(k+1)*Q1(:,1:k);
 
        %Step 7: Compute \til{R}_{k+1}
        Rt2(:,1:k)=Rt2(:,1:k)-conj(alpha(k+1))*Qt1(:,1:k);
        
        k=k+1;
       
        % shift the matrices:
        R1=R2;  Rt1=Rt2;
        Q0=Q1;  Qt0=Qt1; 

        
        result.timing_iteration(k-1)=toc(t0);
    end
    % One extra time Step 8-10, to compute beta(m) and gamma(m)
    omega = this_left_right_scalar_prod(nep,Rt1,R1,m,m);
    beta(m) = sqrt(abs(omega));        
    gamma(m) = conj(omega) / beta(m);

    alpha=alpha(2:end);  % \alpha_1 stored in alpha(2)
    beta=beta(2:end);    % we do not need \beta_1
    gamma=gamma(2:end);  % we do not need \gamma_1

    fprintf('done \n');

    T = spdiags([beta(1:m-1) alpha(1:m-1) [0;gamma(1:m-2)]], -1:1, m-1, m-1);
    result.T=T;
    result.Q_basis = Q_basis;
    result.Qt_basis = Qt_basis;

    
function c=left_right_scalar_prod(nep,At,B,ma,mb)
% Compute the scalar product based on the function nep.M_lin_comb 
% 
    c=0;
    % This is the nasty double loop, which brings
    % complexity O(m^3n). Will be limiting if we do many iterations
   XX=zeros(size(B,1),mb); % pre-allocate
   for j=1:ma
       dd=1./factorial(j:(j+mb-1));
       XX=bsxfun(@times,B(:,1:mb),dd);  % Column scaling: Faster than constructing
                                        % diagonal matrices and multiplying
       z=-nep.M_lin_comb(XX,j);
       c=c+At(:,j)'*z;
   end  
    