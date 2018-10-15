function [V,T,omega] = Lanczos_nep_full(nep,k,q1)
%LANCZOS Indefinite Lanczos method
%   Implementation of the indefinite Lanczos method as 
%   described in the manuscript. Three vertors are used.
%
% TODO:
% Fix now the multiplication with the matrix L

n=nep.n;
% projection matrix
T=zeros(k+1,k);
% initialize vectors of the three term recurrence
Q=zeros(n,k+1); Qp=zeros(n,k+1); Qn=zeros(n,k+1);
Q(:,1)=q1; 
% initialize the vector containing the B-norms
omega=zeros(k+1,1);

omega(1)=Q(:,1)'*(nep.Md(1,Q(:,1)));
% initialize the matrix that will contain the
% first block row of the Krylov basis
V=zeros(n,k+1);
V(:,1)=q1;

% precompute the matrix that factorizes SB
G=zeros(k+1,k+1);
for i=1:k+1 
    G(i,1)=1/i;
end
for j=1:k
    for i=1:k+1
        G(i,j+1)=(G(i,j)*j)/(i+j);
    end
end

Qn=zeros(n,k+1);
Z=zeros(n,k+1);
p=length(nep.f);

for j=1:k
    
    
    % action of A^(-1) B 
    Qn(:,2:j+1)=bsxfun(@rdivide,Q(:,1:j),1:j);
    Qn(:,1)=nep.Mlincomb(ones(j,1),Qn(:,2:j+1));
    Qn(:,1)=-nep.Minv(Qn(:,1));
    
    % B-multiplcation 
    Z=0*Z;
    for t=1:p
        Z(:,1:j+1)=Z(:,1:j+1)+(nep.A{t}*Qn(:,1:j+1))*(G(1:j+1,1:j+1).*nep.FHD{t}(1:j+1,1:j+1));
    end
    
    % orthogonlization (three terms recurrence)
    alpha=sum(sum(bsxfun(@times,conj(Z),Q)));%=Z(:)'*Q(:);
    if j>1
        beta=sum(sum(bsxfun(@times,conj(Z),Qp)));%=Z(:)'*Qp(:);
    end
    gamma=sum(sum(bsxfun(@times,conj(Z),Qn)));%=Z(:)'*W(:);
    
    % from this point nothing changes
    T(j,j)=alpha/omega(j);    
    if j>1
        T(j-1,j)=beta/omega(j-1);
    end
    
    Qn=Qn-T(j,j)*Q;
    if j>1
        Qn=Qn-T(j-1,j)*Qp;
    end
    
    T(j+1,j)=norm(Qn,'fro');% =norm(W_orth(:))
    Qn=Qn/T(j+1,j);
    
    omega(j+1)=gamma-2*T(j,j)*alpha+T(j,j)^2*omega(j);
    if j>1
        omega(j+1)=omega(j+1)-2*T(j-1,j)*beta+T(j-1,j)^2*omega(j-1);
    end
    omega(j+1)=omega(j+1)/T(j+1,j)^2;

    V(:,j+1)=Qn(:,1);       
    % shift the vectors
    Qp=Q;   Q=Qn;
    
end

end