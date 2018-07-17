function [ T, Omega, W, H ] = Infinite_Lanczos_reduced( v1, m, nep )
%Infinite Indefinite Lanczos for nonlinear eigenvalue problems 
%   Benchmark matlab implementation

B_action= @(y)  real(B_action_fun(nep,y));
BB_dot=@(v,w)   real(v'*B_action(w));


n=length(v1);
W=zeros(n,m+1);

vp=zeros(n,1);
vpp=zeros(n,1);


vp(1:n,1)=v1;
vp=vp/norm(vp);
W(:,1)=vp(1:n);

Omega=zeros(m+1);

Omega(1,1)=BB_dot(vp(1:n),vp(1:n));


H=zeros(m+1,m);

% we want to derive a formula for Z (generators of V)

for j=1:m
    j
    jj=1:(j+1)*n;
                       
    % generate next vector (C action)
    x=reshape(vp,n,j);  y=zeros(n,j+1);
    y(:,2:j+1) = bsxfun(@rdivide, x, 1:j);
    y(:,1)=nep.Md_lin_comb(y(:,2:j+1),j);
    y(:,1)=-nep.M_inv(y(:,1));
    zz=reshape(y,n*(j+1),1);
    vn=zz; 
    
    % at each iteration expand one block
    vp=[vp; zeros(n,1)]; vpp=[vpp; zeros(n,1)];
    
    % three term recurrence (orthogonalization))
    z=B_action(vn);
    %z=vn;
    alpha=z'*vp;
    gamma=z'*vn; 
    H(j,j)=alpha/Omega(j,j);    
    vn=vn-H(j,j)*vp;
    Omega(j+1,j+1)=gamma-2*H(j,j)*alpha+H(j,j)^2*Omega(j,j);
    
    if j>1
        beta=z'*vpp;
        H(j-1,j)=beta/Omega(j-1,j-1);
        vn=vn-H(j-1,j)*vpp;
        Omega(j+1,j+1)=Omega(j+1,j+1)-2*H(j-1,j)*beta+H(j-1,j)^2*Omega(j-1,j-1);        
    end

    H(j+1,j)=norm(vn);
    vn=vn/H(j+1,j);
    
    Omega(j+1,j+1)=Omega(j+1,j+1)/H(j+1,j)^2;    


    % shift the vectors
    vpp=vp; vp=vn;
    
    % extend the basis of the projection space
    W(:,j+1)=vn(1:n);
    %W(:,j+1)=W(:,j+1)-W(:,1:j)*(W(:,1:j)'*W(:,j+1));
    %W(:,j+1)=W(:,j+1)-W(:,1:j)*(W(:,1:j)'*W(:,j+1));
    %W(:,j+1)=W(:,j+1)/norm(W(:,j+1));
    
end

T=Omega*H;
[W,~]=qr(W,0);

end


function z=B_action_fun(nep,y)

    n=nep.n;
    m=length(y)/n;
    cc=nep.cc(1:m,1:m);
    [U,S,V] = svd(cc,0); tol=1e-15;  
    s=diag(S);      i=sum(s./s(1)>tol); %i=1;
    V=V(:,1:i);     U=U(:,1:i);     S=S(1:i,1:i);
    U=U*sqrt(S);    V=V*sqrt(S);
    e=ones(n,1);    UU=kron(U,e);   VV=kron(V,e);    
    
    z=sum(BB_action(nep,UU.*(repmat(y,[1 i]))).*VV,2);


end



function [ Z ] = BB_action( nep, Y )
%BB_ACTION Compute the product B*Y

Z=0*Y;
for j=1:size(Y,2)
    Z(:,j)=BB_action_dep( nep, Y(:,j) );
end


end


function [ Y ] = BB_action_dep( nep, y )


n=nep.n;                % size of the problem
m=length(y)/n;          % number of blocks

Y=reshape(y,n,m);
v=zeros(m,1);  for j=1:m; v(j)=(-1)^j; end

vec=@(x) x(:);

yy=Y(:,1);
z=nep.A1*(Y*v);
Y=bsxfun(@times,z,v');
Y(:,1)=Y(:,1)+yy;
Y=vec(Y);
Y=-Y;

end


