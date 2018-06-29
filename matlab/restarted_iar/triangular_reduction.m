function [ P, T ] = triangular_reduction( H, Ritz_vec )
%TRIANGULAR_REDUCTION reduce a matrix H in triangular form
%   Given H a matrix (m+1)xm and Ritz_vec the set of eigenvectors, the
%   function return the matrices T and P such that
%
%   (P' 0) H P = T
%   (0  1)
%
%   The matrix T has the following struncture (almost triangular)
%
%   T =
%   ( s1  *   *   *   *  )  
%   (     s2  *   *   *  )
%   (         s3  *   *  )
%   (             ... *  )
%   (                 sm )
%   ( *   *   *   *   *  )
%
%   where s1,...,sm are the eigenvalues corresponding to
%   Ritz_vec(:,1),...,Ritz_vec(:,m)
%
%   EXAMPLE OF EXECUTION
%   n=5;    H=rand(n+1,n);
%   [Ritz_vec,Ritz_val] = eig(H(1:end-1,1:end));
%   [ P, T ] = triangular_reduction( H, Ritz_vec );
%   PP=eye(n+1);    PP(1:end-1,1:end-1)=P;
%   PP'*H*P
%   error=norm(PP'*H*P-T)
%
%
%   This function is written in a very naive way, it can be done better.
%
%   Giampaolo Mele
%   27/04/2015




Ritz_vec=fliplr(Ritz_vec);
m=size(H,2);
P=eye(m);

T=H;
for k=m:-1:2
    %keyboard
    Ptemp=eye(m);
    Ptemp(m-k+1:end,m-k+1:end)=HOUSEHOLDER(Ritz_vec(m-k+1:end,k));
    Ritz_vec=Ptemp*Ritz_vec;
    PP=eye(m+1);    PP(1:end-1,1:end-1)=Ptemp;
    T=PP*T*Ptemp';
    
    
    P=Ptemp*P;
end

P=P';
T(1:m,1:m)=triu(T(1:m,1:m));

end

function [ P ] = HOUSEHOLDER( y )
%HOUSEHOLDER P*y = alpha e_1
%   this transformation works also for a complex vector y, is the complex
%   version of Householder matrix. For details check it out 
%   http://en.wikipedia.org/wiki/QR_decomposition

n=max(size(y));

e1 = zeros(n,1);
e1(1) = 1;

alpha=-exp(sqrt(-1)*angle(y(1)))*norm(y);
u=y-alpha*e1;
u=u/norm(u);

w=(y'*u)/(u'*y);

P = eye(n) - (1+w)*(u*u');


end

