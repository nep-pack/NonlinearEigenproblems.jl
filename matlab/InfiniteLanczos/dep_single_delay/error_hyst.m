function [ V, D, err ] = error_hyst( Q, T, H, Omega, nep )
%ERROR_HYST Compute error hystory
%   Detailed explanation goes here


m=size(T,2);
n=nep.n;
P=1;
%A0=full(nep.A0); A1=full(nep.A1); II=eye(n);
err=NaN(m,m);
for k=1:P:m
    fprintf('Generation error hystory, step %d of %d \n',k,m);
    
    % compute Ritz pairs
    [Z,D]=eig(T(1:k,1:k),Omega(1:k,1:k));   D=1./diag(D);    V=Q(:,1:k)*Z;
    %[Z,D]=eig(Omega(1:k,1:k),T(1:k,1:k));   D=diag(D);    V=Q(:,1:k)*Z;
    %[Z,D] = eig(H(1:k,1:k)); D=1./diag(D);   V=Q(:,1:k)*Z;
    %HH=balance(H(1:k,1:k)); [~,D] = eig(HH); D=1./diag(D);   

    
    for j=1:k
        err(j,k)=nep.resnorm(D(j),V(1:n,j));
        %err(j,k)=min(svd(-D(j)*II+A0+exp(-D(j))*A1));
    end
    [err(1:k,k), I]=sort(err(1:k,k),'ascend');

    D=D(I);    
    V=V(:,I);    
    
end

figure(1)
mplot='-k';
for i=1:m-1
%    clf
    semilogy(1:P:m-1,err(i,1:P:m-1),mplot)
%    pause(1)
    hold on
end
axis([0 m-1 1e-18 1e0]);

end

