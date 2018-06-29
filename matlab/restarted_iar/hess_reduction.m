function [ Q, H ] = hess_reduction( L )
%HESS_REDUCTION reverse rectangular Hessenberg
%
%   NAIVE IMPLEMENTATION! I AM SURE CAN BE DONE BETTER
%   Giampaolo Mele 24/04/2015
%
%   This function take into account the case where L is triangular (see
%   notes) with few zeros in the last row.
%
%   reverse rectangular Hessenberg.
%   Given a matrix L rectangular matrix L of size (j+1,j), 
%   the output is a square matrix Q such that
%   ( Q' 0 ) * L * Q
%   ( 0  1 ) 
%   is an Hessenberg matrix with positive elements in the subdiagonal. 
%   Everithing is done with Householder matrices.
%   This algorithm works also for complex matrices, for details see my 
%   master thesis 
%   or 
%   http://web.eecs.utk.edu/~dongarra/etemplates/node295.html
%   or 
%   A. Ruhe. Eigenvalue algorithms with several factorizations - 
%   a unified theory yet? Technical Report 1998:11, 
%   Department of Mathematics, Chalmers University of Technology, 
%   Göteborg, Sweden, 1998. 
%
%   EXAMPLE OF EXECUTION
%
%
%   n=50;
%   L=rand(n+1,n)+rand(n+1,n)*1i;
%   [ Q ] = hess_reduction( L );
%   QQ=eye(n+1);
%   QQ(1:n,1:n)=Q;
%   H=QQ'*L*Q;
%   spy(H)
%   error=norm(QQ*H*Q'-L)
%
% 
%   n=100;
%   L=rand(n+1,n)+rand(n+1,n)*1i;
%   L=triu(L);
%   L(end,:)=rand(1,n);
%   L(end,1:10)=0;
%   [ Q, H ] = hess_reduction( L );
%   QQ=eye(n+1);
%   QQ(1:n,1:n)=Q;
%   spy(H)
%   error=norm(QQ*H*Q'-L)
%
%
%   Note on particular structures: (in Arnoldi method appears)


n=min(size(L));

zeros_count=0;
while(L(end,zeros_count+1)==0)
    zeros_count=zeros_count+1;
end

if(zeros_count==0)

    T1 = fliplr(eye(n+1));
    T2 = fliplr(eye(n));

    P = TRIANG( (T1*L*T2)' );

    P=T2*P*T2;

    PP=eye(n+1);
    PP(1:n,1:n)=P;

    [ Q ] = POSITIVE_HESS( PP'*L*P );

    Q=P*Q;


    QQ=eye(n+1);
    QQ(1:n,1:n)=Q;
    H=QQ'*L*Q;
    H=triu(H,-1);   % enforce Hessenberg structure

else
    [ Qtemp, ~ ] = hess_reduction( L(zeros_count+1:end,zeros_count+1:end) );
    Q=eye(n);
    Q(zeros_count+1:end,zeros_count+1:end)=Qtemp;
    QQ=eye(n+1);    QQ(1:n,1:n)=Q;
    
    H=QQ'*L*Q;
    H=triu(H,-1);   % enforce Hessenberg structure
    
end
    



end

function [ P ] = TRIANG( A )
% TRIANG triangular form.
%   given a rectangular matrix A of size (j,j+1) the output is a matrix P 
%   such that 
%   P' * A * ( 1 0 )
%            ( 0 P )
%   is a triangular matrix

j = min(size(A));

P = eye(j);

for i = 1:1:j-1
    
    Q=eye(j);
    Q(i:j,i:j) = HOUSEHOLDER( A(i:j,i) );
    P=Q*P;
    Q2=eye(j+1);
    Q2(2:j+1,2:j+1)=Q';
    A=Q*A*Q2;

end

P=P';

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

function [ P ] = POSITIVE_HESS( H )
%POSITIVE_HESS Positive Hessenberg form
%   Given an Hessenberg matrix H with complex numbers in the subdiagonal
%   this function return a phase matrix P (diagonal matrix with 
%   unitary elements) such that
%
%   ( P' 0 ) H P
%   ( 0  1 )
%

%tol=1.e-14;
tol=eps;

j=size(H, 2);

P=eye(j);


for i=j:-1:1
    
    Q=eye(j);
    Q2=eye(j+1);

    
    if(abs(imag(H(i+1,i)))>tol)
        Q(i,i) = exp(-angle(H(i+1,i))*sqrt(-1));
        Q2(1:j,1:j) = Q;
        H=Q2'*H*Q;
        P=P*Q;
    end
    
    if(real(H(i+1,i))<0)
        Q(i,i) = -Q(i,i);
        Q2(1:j,1:j) = Q;
        H=Q2'*H*Q;
        P=P*Q;
    end
    
    
    
    

end

end
