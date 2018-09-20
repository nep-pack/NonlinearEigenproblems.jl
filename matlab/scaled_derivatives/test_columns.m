close all
clear all
clc
%% DEFINE THE PROBLEM
n=5;                % problems size (wanted n-1 scaled-derivative)
alpha=rand(1,n);    % rescaling factors
mu=rand;            % the derivatives are evaluated in mu

%% COMPUTE THE SCALED DERIVATIVES WITH A MATRIX FUNCTION EVALUATION
JJ=mu*diag(ones(n,1))+diag((1:1:n-1).*alpha(2:n)./alpha(1:n-1),-1);  % tilde J

e1=eye(n,1);
FF=alpha(1)*funm(JJ,@myfun)*e1;


%% COMPUTE THE SCALED-DERIVATIVES IN THE STANDARD WAY
F=zeros(n,1);
for k=1:n
    F(k)=myfun(mu,k-1)*alpha(k);
end


%% 	COMPARE THE RESULTS
F
FF
norm(F-FF)