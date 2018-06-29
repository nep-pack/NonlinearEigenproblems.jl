function [ P ] = Pinv_mat( n, rho, gamma )
%P Summary of this function goes here
%   Detailed explanation goes here
I=eye(n);
P=zeros(n);
for j=1:n
    P(:,j)=mon2cheb(rho,gamma,I(:,j));
end
end