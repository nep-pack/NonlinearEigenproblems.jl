function [ V, H ] = InfArn( nep, v1, m, variation )
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




% iterations of the algorithm
for k=1:m
    fprintf("Iteration %d\n",k);
    % COMPUTING NEXT VECTOR OF THE ARNOLDI SEQUENCE
    
    y=zeros(n,k+1);   % vector-coefficients that defines the next vector
                      % of the Arnoldi sequence
    
    % reshape the last vector of the Arnoldi sequence 
    % (see article about waveguides eigenproblem)
    
    % computing y2,...,y_{k+1}  
    x=reshape(V(:,k),n,k);
    if strcmp(variation,'Taylor')
        for j=2:k+1        
           y(:,j)=1/(j-1)*x(:,j-1);        
        end
    elseif strcmp(variation,'Chebyshev')
        y(:,2:end)=x*LL(k);
    else
        print('Error')
        break
    end
    
    % computing y1  
    if strcmp(variation,'Taylor')
        y(:,1)=zeros(n,1);    
        for s=1:k
            y(:,1) = y(:,1) + nep.Mdd(s)*y(:,s+1);
        end
        y(:,1)=-M0\y(:,1);
    elseif strcmp(variation,'Chebyshev')
          y(:,1)=compute_y0(x,y,nep);
    else
        print('Error')
        break          
    end    
    y=reshape(y,(k+1)*n,1);
    
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
function y0=compute_y0(x,y,nep)
    tt=1; % delay
    n=length(nep.A0);
    N=size(x,2);
    k=2/(nep.b-nep.a) ; c=(nep.a+nep.b)/(nep.a-nep.b);

    
    % Chebyshev polynomials of the first kind
    T=@(n,x) cos(n*acos(x));
    TT=@(n,x) T(n,k*x+c);
    
    % Chebyshev polynomials of the second kind
    %U=@(n,x) sin((n+1)*acos(x))/sin(acos(x));
    U=@(n,x) n+1; % observe that U(n,1)=n+1 and we evaluate U only in 1 in this code
    
    y0=zeros(n,1);
    for i=1:N+1
        y0=y0-y(:,i);
    end
    
    for i=1:N
        y0=y0+nep.A0*x(:,i);
    end
    
    for i=1:N-1
        y0=y0+nep.A1*T(i,1-2*tt)*x(:,i+1);
    end
    y0=(nep.A0+nep.A1)\y0;
    
end


