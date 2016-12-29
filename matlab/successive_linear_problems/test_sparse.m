close all
clear all
clc

n=10;
A0=rand(n); A1=rand(n); I=eye(n);

M =@(s) s*I+A0+A1*exp(-s);
Mp=@(s) I-A1*exp(-s);


s=0;
err=[];
for i=1:10
    [V,D] = eigs(M(s),Mp(s),1,s);
    D=diag(D);
    [~,jj] = min(abs(D-s));
    s=s-D(jj);
    v=V(:,jj);
    err=[err norm(M(s)*v)];
end

semilogy(err,'--*')