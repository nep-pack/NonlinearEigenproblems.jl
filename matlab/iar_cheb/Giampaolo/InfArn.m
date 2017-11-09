function [ V, H ] = InfArn( nep, v1, m )
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



% iterations of the algorithm
for k=1:m
    
    % COMPUTING NEXT VECTOR OF THE ARNOLDI SEQUENCE
    
    y=zeros(n,k+1);   % vector-coefficients that defines the next vector
                      % of the Arnoldi sequence
    
    % reshape the last vector of the Arnoldi sequence 
    % (see article about waveguides eigenproblem)
    
    Vsplit=reshape(V(:,k),n,k);
                      
    for j=2:k+1        
        y(:,j)=1/(j-1)*Vsplit(:,j-1);        
    end
    
    % computing y1
    
    y(:,1)=zeros(n,1);    
    for s=1:k
        y(:,1) = y(:,1) + nep.Mdd(s)*y(:,s+1);
    end    
    y(:,1)=-M0\y(:,1);
    y=reshape(y,(k+1)*n,1);
    
    % expand V
    V=[V ; zeros(n,k)];
    
    % orthogonalization
    H(:,k)=V'*y;
    y=y-V*H(:,k);
    
    
    H(k+1,k)=norm(y);
    V(:,k+1)=y/H(k+1,k);
    
    
    
    
end



    
    
    
end