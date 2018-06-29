function [ V, H ] = infarn_exp(nep, C, S, Y, kmax, pl )
%INFARN_EXP Infinite Arnoldi algorithm exponential version
%   Perform kmax+1 steps of infinite Arnoldi algorithm starting from an
%   exponential function, that is the first starting function is
%   phi(theta)=Y*exp(S*theta)c and a locked part of size pl
%
%   Giampaolo Mele
%   26/11/2015
%
%   Additional notes:
%   Reorthogonalization: YES
%   Preallocation of variables: YES
%   Compression of V: NO
%   LU-factorization of M(0): YES
%   Optimizations: NONE

   
n=size(Y,1);                % size of the problem
V=zeros(n*kmax,kmax+1);     % polynomial part
H=zeros(kmax,kmax-1);       % initizization of H
p=size(Y,2);                % restarting paramenter
tol=eps;                    % tollerance (converging series)

% exstracting the converged part
if pl~=0
	H(1:pl,1:pl)=S(1:pl,1:pl)^(-1);
end

% normalization of the starting function
beta=norm_exp(C(:,pl+1),S,Y); 
C(:,pl+1)=C(:,pl+1)/beta;

for k=pl+1:kmax
      
    
    % EXPONENTIAL PART
    C(:,k+1)=S\C(:,k); 
    
    
     % COMPUTE x
     v=reshape(V(1:(k-pl-1)*n,k),n,k-pl-1);
     x=zeros(n,k-pl);
     
     for i=1:k-pl-1
         x(:,i+1)=v(:,i)/i;    
     end
     
     for i=1:k-pl-1
         x(:,1)=x(:,1)+T(i,nep)*x(:,i+1);
     end
   
     res=1;
     i=k-pl;
     while norm(res)>tol
         res=T(i, nep)*Y*S^(i)*C(:,k+1)/factorial(i);
         x(:,1)=x(:,1)+res;
         i=i+1;
     end
          
     x(:,1)=-nep.Q*(nep.U\(nep.L\(nep.P*((x(:,1))./nep.R))));
        
     x=x(:);
    
     % EXPANSION OF THE POLYNOMIAL PART
     V((k-pl-1)*n+1:((k-pl-1)+1)*n,1:k)=Y*S^(k-pl-1)*C(:,1:k)/factorial(k-pl-1);
    
     % ORTHOGONALIZATION (compute w_orth)
     % polynomial part
     hp=V(1:(k-pl)*n,1:k)'*x;
   
     % exponential part
     res=1;    sum=zeros(p,p);
    
     i=k-pl;    
     while norm(res)>tol
         res=Y*S^(i);
         res=res'*res;
         res=res/(factorial(i)^2);        
         sum=sum+res;
         i=i+1;
     end
     he=C(:,1:k)'*sum*C(:,k+1);
    
     h=hp+he;
    
     x_orth=x-V(1:(k-pl)*n,1:k)*h;
     C(:,k+1)=C(:,k+1)-C(:,1:k)*h;
    
     % REORTHOGONALIZATION
     % polynomial part
     hhp=V(1:(k-pl)*n,1:k)'*x_orth;
   
     % exponential part
     res=1;    sum=zeros(p,p);
    
     i=k-pl;    
     while norm(res)>tol
         res=Y*S^(i);
         res=res'*res;
         res=res/(factorial(i)^2);        
         sum=sum+res;
         i=i+1;
     end
     hhe=C(:,1:k)'*sum*C(:,k+1);
    
     hh=hhp+hhe;    
    
     x_orth=x_orth-V(1:(k-pl)*n,1:k)*hh;
     C(:,k+1)=C(:,k+1)-C(:,1:k)*hh;    
            
     % CONTRIBUTION OF THE REORTHOGONALIZATION
     h=h+hh;
    
     % COMPUTE THE NORM
     beta_p=x_orth'*x_orth; % polynomial
     
     res=1;
     i=k-pl;
     
     beta_e=0;  % exponential
     while norm(res)>tol
         res=Y*S^(i)*C(:,k+1);  
         res=res'*res;          
         res=res/(factorial(i)^2);
         beta_e=beta_e+res;        
         i=i+1;
     end    
     beta=sqrt(beta_p+beta_e);
    
    
     % NORMALIZATION
     V(1:(k-pl)*n,k+1)=x_orth/beta;
     C(:,k+1)=C(:,k+1)/beta;
    
     % UPDATE THE MATRIX H
     H(1:k,k)=h; H(k+1,k)=beta;
    
    
    
end

V=V(1:n,:);


end








function beta=norm_exp(c,S,Y)
% norm_exp compute the norm of the sum of an exponential function 
% of the form phi(theta)=Y*exp(S*theta)c


    % tolerance for the truncation of the exponential part
    tol=eps;
    
    beta=0;
 
    
    res=1;
    i=0;
    
    while norm(res)>tol
        
        res=Y*S^(i)*c;  
        res=res'*res;  
        
        res=res/(factorial(i)^2);
        beta=beta+res;
        
        i=i+1;
    end
    
    beta=sqrt(beta);
    

    
end




function [ val ] = T( k, nep )
%T Derivatives in zero of the DEP
%   The DEP (Delay Eigenvalue Problem) is defined by
%   T(lambda) = -lambda*I+A_1*exp(-lambda)+A_2


    n=length(nep.A1);
    
    
    if (k==0)
        val=nep.A1+nep.A2;
    elseif(k==1)
        val=-eye(n)-nep.A1;
    else
        val=nep.A1*(-1)^k;
    end


end
