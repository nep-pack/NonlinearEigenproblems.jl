function [ T, Omega, W, H ] = Infinite_Lanczos( v1, m, nep )
%Infinite Indefinite Lanczos for nonlinear eigenvalue problems 
%   Benchmark matlab implementation


B_action= @(y)  B_action_fun(nep,y);



BB_dot=@(v,w)   real(v'*B_action(w));

BB_dott=@(v,w) real(B_dott( nep, w, v ));  % (indefinite) scalar product


n=length(v1);
W=zeros(n,m+1);

vp=zeros((m+1)*n,1);
vpp=zeros((m+1)*n,1);
vn=zeros((m+1)*n,1);

V=zeros((m+1)*n,3);


V(1:n,2)=v1;
V(:,2)=V(:,2)/norm(V(:,2));
W(:,1)=V(1:n,2);

Omega=zeros(m+1);

Omega(1,1)=BB_dot(vp(1:n),vp(1:n));

H=zeros(m+1,m);

for j=1:m
                       
    j
    % generate next vector
    %V(1:(j+1)*n,3)=C_action(nep,V(1:j*n,2));
    
    Y=reshape(V(1:j*n,2),n,j);
    Y(:,2:j+1) = bsxfun (@cdivide, Y(:,1:j), 1:j);    
    Y(:,1)=Md_lin_comb(Y,j);
    Y(:,1)=-M_inv(Y(:,1));
    vvn(1:(j+1)*n)=Y(:);
    
    norm(vn-vvn)
    pause
    
    jj=1:(j+1)*n;
    
    % B-orthogonalize
    h=zeros(j,1);
    %norm(BB_dot(vn(jj),vp(jj))-BB_dott(vn,vp))   
    h(j)=Omega(j,j)\BB_dot(vp(jj),vn(jj));
    vn(1:(j+1)*n)=vn(1:(j+1)*n)-vp(1:(j+1)*n)*h(j);        
    
    
    if j>1
        %norm(BB_dot(vn(jj),vpp(jj))-BB_dott(vn,vpp))           
        h(j-1)=Omega(j-1,j-1)\BB_dot(vpp(jj),vn(jj));
        vn(1:(j+1)*n)=vn(1:(j+1)*n)-vpp(1:(j+1)*n)*h(j-1);   
    end
    
    % B-re-orthogonalize
    % MAYBE NOT NEEDED
    
    H(1:j,j)=h;
    
    H(j+1,j)=norm(vn);
    vn=vn/H(j+1,j);

    %norm(BB_dot(vn(jj),vn(jj))-BB_dott(vn,vn))   
    %pause
    
    Omega(j+1,j+1)=BB_dot(vn(jj),vn(jj));
    

    % shift the vectors
    vpp=vp;
    vp=vn;
    
    % extend the basis (maybe orth at each iter, fix later)
    W(:,j+1)=vn(1:n);
    W(:,j+1)=W(:,j+1)-W(:,1:j)*(W(:,1:j)'*W(:,j+1));
    W(:,j+1)=W(:,j+1)-W(:,1:j)*(W(:,1:j)'*W(:,j+1));
    W(:,j+1)=W(:,j+1)/norm(W(:,j+1));
    
    %norm(W(:,1:j+1)'*W(:,1:j+1)-eye(j+1))
    
    
    
end

T=Omega*H;


end


function z=B_action_fun(nep,y)

    n=nep.n;
    m=length(y)/n;
    cc=nep.cc(1:m,1:m);
    [U,S,V] = svd(cc,0); tol=1e-15;  
    s=diag(S);      i=sum(s./s(1)>tol); 
    %i=min(6,i);
    V=V(:,1:i);     U=U(:,1:i);     S=S(1:i,1:i);
    U=U*sqrt(S);    V=V*sqrt(S);
    e=ones(n,1);    UU=kron(U,e);   VV=kron(V,e);
    %B_action= @(y)  sum(BB_action(nep,UU.*y).*VV,2);
    z=(BB_action2(nep,UU.*(y*ones(1,i))).*VV)*ones(i,1);


end


function [ y ] = C_action( nep, x )
%C_action compute C*x where C is the companion linearization
%   compute the action of the companion matrix C in a block-structured
%   vector x

n=nep.n;                % size of the problem
k=length(x)/n;          % number of blocks


x=reshape(x,n,k);
y=zeros(n,k+1);


for j=2:k+1        
    y(:,j)=1/(j-1)*x(:,j-1);        
end

for s=1:k
    y(:,1) = y(:,1) + nep.Md(s)*y(:,s+1);
end

y(:,1)=-nep.M_inv(y(:,1));
y=y(:);


end

function [ gamma ] = B_dott( nep, x, y )
%B_ACTION Compute the product x'*B*y 
%   The complexity of this function depends of the number of factors p
%   needed to write the nonlinear eigenvalue problem with
%   M(l) = \sum_{s=1}^p f_s(l) A_s
%   The larger the p, the slower the method


n=nep.n;                % size of the problem
m=length(y)/n;          % number of blocks

y=reshape(y,n,m);
x=reshape(x,n,m);


% precomputation
xA1y=x'*(nep.A1*y);
xy=x'*y;

xMdy1=-xy-xA1y;
xMdy=xA1y;

g=zeros(2*m+1,1);
for k=2:2*m+1
    g(k)=(-1)^k;
end
             
gamma=0;
for i=1:m
    for j=1:m 
        gamma=gamma+nep.cc(i,j)*(xMdy(i,j)*g(j+i-1)+xMdy1(i,j)*(j+i-1==1));
    end
end

end
