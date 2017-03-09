randn('seed',0)
A=randn(100);  

d=eig(A);
A=A-eye(length(A))*(d(1)+eps*100);

%delta=sqrt(eps);
%A=A-delta*speye(size(A));
b=randn(100,1);
n=100;
kmax=70;
Q=zeros(n,kmax);
H=zeros(kmax+1,kmax);
Q(:,1)=b/norm(b);
nv=[];
for k=1:kmax
    w=A*Q(:,k);

    Q1=Q(:,1:k);
    h=Q1'*w;
    w2=w-Q1*h;
    g=Q1'*w2;
    w2=w2-Q1*g;
    h=g+h;

    beta=norm(w2);
    h2=[h;beta];
    H(1:(k+1),k)=h2;
    Q(:,k+1)=w2/beta;
    
    % Extract by svd of Hessenberg matrix
    z=H(1:k+1,1:k)\eye(k+1,1);
    
    %[U,S,V]=svd(H(1:k,1:k));
    %z=V(:,end);
    v=Q(:,1:k)*z; 
    v=v/norm(v);
    nv=[nv,norm(A*v)];
end

semilogy(nv,'r')
%z=H\eye(k+1,1);
%v=Q(:,1:k)*z; 

