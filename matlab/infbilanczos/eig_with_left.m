function [V,D,W]=eig_with_left(A)
% Computes left and right eigenvectors 
%
%   A=randn(5);
%   [V,D,W]=eig_with_left(A);
%   norm(W'*A-D*W')
%
%   same as [V,D,W]=eig(A) in newer versions of matlab

    [V,D]=eig(A);
    d=diag(D);
    [WW,DD]=eig(A');
    dd=diag(DD);
     
    W=zeros(size(D));
    for i=1:size(WW,2)
        [Y,I]=min(abs(d(i)-dd')); 
        W(:,i)=WW(:,I);
    end