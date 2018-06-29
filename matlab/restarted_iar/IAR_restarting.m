function [ ERR_restarting ] = IAR_restarting( nep, S, Y, kmax, p, N_restart )
%IAR_restarting Infinite Arnoldi with Restart for a DEP
%
%   INPUT
%
%   The matrices nep.A1 and nep.A2 defines the Delay Eigenvalue Problem
%    -lambda*x + nep.A2*x + exp(-lambda)*nep.A1*x = 0
%   
%   m is the number of iterations of IAR before the restart
%
%   p is the number of Ritz Pairs that will be taken in the restart
%   if a Ritz Pair converged, then will be locked after the restart
%   
%   S is an estimation of one eigenvalue (at least can be choosen randomly
%   in the region of interest)
%
%   Y is an estimation of the eigenvector corresponding to S. In particular
%   the eigenvector of the infinite dimensional linear eigenvalue problem
%   will be exp(SY), see notes.
%
%   N_restart is the number of restart that will be performed
%
%   OUTPUT
%
%   ERR_restarting contains the error history


% inizialization (start with an exponential function)
C=1;

% is the parameter used to lock converged Ritz Pairs. That is if the
% residual is less than ltol the Ritz Pair will be locked.
ltol=1.e-14;

% initialization Error History
ERR=NaN(kmax,1);   ERR_restarting=NaN(kmax,kmax,N_restart);

% inizialization converged Ritz Pairs
pl=0;   % zero Ritz Pairs converged

% size of the problem
n=length(nep.A1);   

for restart=1:N_restart
    
    % run the Arnoldi algorithm with m steps
    [ V, H ] = infarn_exp(nep, C, S, Y, kmax, pl );   
    
    
    % compute Ritz pairs
   [Ritz_vec, Ritz_val] = eig(H(1:end-1,1:end));    Ritz_val=diag(Ritz_val);

    % compute eigenvectors approximations
    Eig_vec_approx=V(:,1:kmax)*Ritz_vec;   Eig_vec_approx=Eig_vec_approx(1:n,:);
     
     % sort the approximations of the eigenvectors (according to the residual)
     for i=1:kmax
         Eig_vec_approx(:,i)=Eig_vec_approx(:,i)/norm(Eig_vec_approx(:,i));
         ERR(i)=backward_error(nep, 1/Ritz_val(i), Eig_vec_approx(:,i) );
     end
     [ERR, I]=sort(ERR);
     Ritz_vec=Ritz_vec(:,I);
     
     
     % computing the error history
     ERR_restarting(:,:,restart)=error_history_restarting_version(nep, V, H);
     
     % compute the triagular reduction (see notes)
     [ P, T ] = triangular_reduction( H, Ritz_vec ); 
     PP=eye(length(P)+1);    PP(1:end-1,1:end-1)=P;
     
     % computing the locked part
     pl=0;   % locked part
     while(abs(T(end,pl+1))<ltol)
             T(end,pl+1)=0;     % lock the eigenvalue
             pl=pl+1;
     end
     
     if pl>=p
         error('Increase p')
     end
     
     % first change of basis
     V=V*PP;
     
     % redoucing the number of vectors in the Arnoldi sequence
     V=[V(:,1:p) V(:,end)];
     T=[T(1:p,1:p); T(end,1:p)];
     
     % Hessenberg reduction (see notes)
     [ P, H ] = hess_reduction( T );
     PP=eye(length(P)+1);    PP(1:end-1,1:end-1)=P;
     
          
     % second change of basis
     V=V*PP;
     
     
     % new guess (invariant pair approximation)
     S=H(1:p,1:p)^(-1); Y=V(1:n,1:p);   

     C=eye(p,kmax);    C(:,pl+2:end)=zeros(size(C(:,pl+2:end)));


     % REORTHOGONALIZATION LOCKED PART
     for kk=1:pl
         [C(:,1:kk),h,~]=GS_exp(Y,S,C(:,1:kk));
         [C(:,1:kk),g,beta]=GS_exp(Y,S,C(:,1:kk));
         h=h+g; C(:,kk)=C(:,kk)/beta;
         H(1:kk-1,kk)=H(1:kk-1,kk)+h;
     end
     

     % if at least one Ritz pair converged, after the lock choose the new
     % function (see notes) and double orthogonalize
     if(pl>0)
         [C(:,1:pl+1),~,~]=GS_exp(Y,S,C(:,1:pl+1));  
         [C(:,1:pl+1),~,beta]=GS_exp(Y,S,C(:,1:pl+1));  
         C(:,pl+1)=C(:,pl+1)/beta;       
     end
     
     
     
     
     
end     
     


end


function ERR_hist=error_history_restarting_version(nep, V, H)


m=size(H,2);
n=size(nep.A1,1);

ERR_hist=NaN(m,m);
for k=1:m
    
    clc
    fprintf('Computation error history \n');
    fprintf('PROGRESS %d %% \n',round((k*100)/m));
    
    [RitzVec, RitzVal] = eig(H(1:k,1:k));
    RitzVal=diag(RitzVal); RitzVal=1./RitzVal;
    RitzVec=V(:,1:k)*RitzVec;   RitzVec=RitzVec(1:n,:);
    
    
    % COMPUTATION OF THE ERROR
    for i=1:k
        RitzVec(:,i)=RitzVec(:,i)/norm(RitzVec(:,i));
        ERR_hist(i,k)=backward_error(nep, RitzVal(i), RitzVec(:,i) );
    end
    
    % SORT ERROR VECTOR
    ERR_hist(1:k,k)=sort(ERR_hist(1:k,k));


end

end


function [ err ] = backward_error(nep, lambda, x )
%backward_error compute the backward error of the following DEP
%    -lambda*x + nep.A2*x + exp(-lambda)*nep.A1*x = 0

    
    err=norm( -lambda*x + nep.A2*x + exp(-lambda)*nep.A1*x);
    err=err/(abs(lambda) + norm(nep.A2,1) + abs(exp(-lambda))*norm(nep.A1,1));
    

end






function beta=norm_exp(c,S,Y)
% norm_exp compute the norm of an exponential function 
%
%   INPUT
%
%   d is the degree of the polynomial part
%
%   x, c, S, Y represent the function/vector which we want to compute the norm
%   x represent the polynomial part   
%   c, S, Y represent the exponential part
%
%   OUTPUT
%
%   beta is the norm of the vector




    % tolerance for the truncation of the exponential part
    tol=eps;
    
    
    res=1;
    i=0;
    beta=0;
    while norm(res)>tol
        
        res=Y*S^(i)*c;  
        res=res'*res;  
        
        res=res/(factorial(i)^2);
        beta=beta+res;
        
        i=i+1;
    end
    
    beta=sqrt(beta);
    

    
end


function [C,h,beta]=GS_exp(Y, S, C)


    % tolerance for the truncation of the exponential part
    tol=1.e-16;
    
    % exponential part
    res=1;
    
    p=size(Y,2);
    sum=zeros(p,p);
    

    i=0;    
    while norm(res)>tol
        res=Y*S^(i);
        res=res'*res;
        res=res/(factorial(i)^2);
        
        sum=sum+res;
        i=i+1;
    end
    
    h=C(:,1:end-1)'*sum*C(:,end);
    
    % orthogonalization
    C(:,end)=C(:,end)-C(:,1:end-1)*h;
    
    
    % computation of the norm
    i=0;
    res=1;
    beta=0;
    while norm(res)>tol
        
        res=Y*S^(i)*C(:,end);  
        res=res'*res;  
        
        res=res/(factorial(i)^2);
        beta=beta+res;
        
        i=i+1;
    end
    
    beta=sqrt(beta);
    
end
