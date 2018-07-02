function [ V, H ] = InfArn_change_basis( nep, v1, m )
%INFARN naive implementation
%   Naive imprementation of Infinite Arnoldi
%   Date: 13 May 2014
%   Giampaolo Mele
%   
%   INPUT
%   M is the matrix-function that defines the NLEP
%   m is the number of steps of the algorithm
%
%   OUTPUT
%   Ritz is the vector containig the Ritz values

% inizialization
M0 = nep.Mdd(0);
n=size(M0,1);   % size of the problem

V=v1/norm(v1);                  % basis of Krylov space
H=zeros(1,1);                   % Hessenberg matrix

LL = @(k) L(k)*(nep.b-nep.a)/4;
nep.LL=LL;

a=nep.a; b=nep.b;

P=P_mat(m+1,2/(b-a),(a+b)/(a-b));
Pinv=Pinv_mat(m+1,2/(b-a),(a+b)/(a-b));

% iterations of the algorithm
for k=1:m
    fprintf("Iteration %d\n",k);
    % COMPUTING NEXT VECTOR OF THE ARNOLDI SEQUENCE
    

    % pep version (to do)
    yy=zeros(n,k+1);
    xx=reshape(V(:,k),n,k);
    yy(:,2:end)=xx*LL(k);
    
    yy(:,1)=pep_y0(xx,yy,nep);
        
    % --------------------------------------------
    % new variation
    y=zeros(n,k+1);
    x=reshape(V(:,k),n,k);
    % change basis (to Monomial)
    x=x*P(1:k,1:k)';
    
    % apply the operator B (Taylor)
    for j=2:k+1        
        y(:,j)=1/(j-1)*x(:,j-1);
    end

    y(:,1)=zeros(n,1);    
    for s=1:k
        y(:,1)=y(:,1)+nep.Mdd(s)*y(:,s+1);
    end
    y(:,1)=-M0\y(:,1);   
    % change basis (to Chebyshev)
    y=y*Pinv(1:k+1,1:k+1)';    
    y=reshape(y,(k+1)*n,1);
    % --------------------------------------------
    
   
    
    % COMPARISON
    yy=reshape(yy,(k+1)*n,1);
    [y(1:n) yy(1:n) y(1:n)-yy(1:n)]
    y0diff=norm(y(1:end)-yy(1:end))
    ynndiff=norm(y(n+1:end)-yy(n+1:end))
    pause
    
    
    
    % expand V
    V=[V ; zeros(n,k)];
    
    % double orthogonalization
    h=V'*y;    y=y-V*h;
    g=V'*y;    y=y-V*g;
    H(:,k)=h+g;
    
    
    H(k+1,k)=norm(y);
    V(:,k+1)=y/H(k+1,k);
    
    
    
    
end



end

% matrix needed to the expansion
function Lmat=L(k)
    if k==1
        Lmat=2;
    else
        Lmat=diag([2, 1./(2:k)])+diag(-1./(1:(k-2)),-2);
    end
end


function y0=pep_y0(x,y,nep)
    a=nep.a;    b=nep.b;
    c=(a+b)/(a-b);
    N=size(x,2);
    Tc=cos((0:N)*acos(c)).';
    n=nep.n; 
    d=length(nep.coeff)-1;
    % sum over all the matrix coefficients
    y0=zeros(n,1);
    
    % compute the derivation matrix
    R=[zeros(1,N); inv(nep.LL(N))]; R=R(1:N,1:N);
    % sum for every coefficient
    v=Tc(1:N);
    for j=0:d-1
        Bj=-nep.coeff{j+2};
        y0=y0+Bj*(x*v);
        v=R*v;
    end
    y0=nep.M0solver(y0);
    
    y0=y0-y*Tc;
    
    
    
    
    
end