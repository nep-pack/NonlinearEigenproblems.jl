function [ P ] = P_mat( n, rho, gamma )
%P Summary of this function goes here
%   Detailed explanation goes here
I=eye(n);
P=zeros(n);
for j=1:n
    P(:,j)=cheb2mon(rho,gamma,I(:,j));
end
end