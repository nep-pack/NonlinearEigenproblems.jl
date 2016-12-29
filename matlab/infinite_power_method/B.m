function [ y ] = B( nep, x )
%B Summary of this function goes here
%   Detailed explanation goes here

% number of the block rows
n=nep.n;    k=length(x)/n;

x=reshape(x,n,k);

y=zeros(n,k+1);

for i=1:k
    y(:,1)=y(:,1)+nep.Md(i,x(:,i)/i);
end
y(:,1)=-nep.M0solver(y(:,1));


for j=1:k
    y(:,j+1) = x(:,j)/j;
end

y=y(:);

end