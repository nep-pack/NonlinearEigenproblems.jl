close all
clear all
clc

% define a polynomial in monomial basis with random coefficients "a"
n=5;    % degree+1
a=rand(n,1);

% % coefficients that one may expect from IAR
%a(:,1)=exp(-(1:1:n));

% set scaling and shift
rho=2;    gamma=3;

% convert the coefficients to Chebishev basis
c=mon2cheb(rho,gamma,a);

% convert the coefficients to Chebishev basis with the naive approach
cc=naive_mon2cheb(rho,gamma,a);

fprintf("Error in comparison to the naive approach %e\n",norm(c-cc));

% re-convert the coefficients to the monomial basis
aa=cheb2mon(rho,gamma,c);

fprintf("Error before and after a loop conversion %e\n",norm(a-aa));

% % UNCOMMENT TO PLOT THE POLYNOMIALS IN BOTH BASIS (SHOULD BE IDENTICAL)
% x=linspace(-1,1,100);
% yC=cheb_polyval(x,c,rho,gamma);
% plot(x,yC,'--r')
% 
% hold on
% y=polyval(flip(a),x);
% plot(x,y,'-k')
% 
% figure
% semilogy(x,abs(y-yC))