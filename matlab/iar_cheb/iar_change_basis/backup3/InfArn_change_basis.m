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
a=nep.a; b=nep.b;



% iterations of the algorithm
for k=1:m
    fprintf("Iteration %d\n",k);
    % COMPUTING NEXT VECTOR OF THE ARNOLDI SEQUENCE
    


    
    % --------------------------------------------
    % new variation
    y=zeros(n,k+1);
    x=reshape(V(:,k),n,k);
    % change basis (to Monomial)
    for i=1:size(x,1)
        x(i,:)=cheb2mon(2/(b-a),(a+b)/(a-b),x(i,:));
    end
    
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
    for i=1:size(y,1)
        y(i,:)=mon2cheb(2/(b-a),(a+b)/(a-b),y(i,:));
    end
    y=reshape(y,(k+1)*n,1);
    % --------------------------------------------

    % --------------------------------------------
    % standard
    yy=zeros(n,k+1);
    xx=reshape(V(:,k),n,k);
    yy(:,2:end)=xx*LL(k);
    yy(:,1)=compute_y0(xx,yy,nep);
    yy=reshape(yy,(k+1)*n,1);
    % --------------------------------------------

    
    
    norm(y-yy)
    
    
    

    
    
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

% function for compuing y0 for this specific DEP
function y0=compute_y0(X,Y,nep)

a=nep.a;    b=nep.b;
A1=nep.A1;  A0=nep.A0;
tau=1;

k=2/(b-a) ; c=(a+b)/(a-b);
cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0= (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));
    
end


% function for compuing y0 for this specific DEP (Elias)
% function y0=compute_y0(X,Y,nep)
% 
% a=nep.a;    b=nep.b;
% A1=nep.A1;  A0=nep.A0;
% tau=1;
% 
% k=2/(b-a) ; c=(a+b)/(a-b);
% cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
% cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
% Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));
% 
% % QDEP: 
% y0= (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));
%     
% end


