close all
clear all
clc

m=5;
C=gen_matrix_c( 2*m );
G=zeros(2*m,2*m);
G(1,:)=1./(1:(2*m));  G(:,1)=1./(1:(2*m));  
for j=1:2*m   % later try to avoid for loops
    for i=2:m
        G(i,j)=C(i-1,j)/j;
    end
end
G=G(1:m,1:m);
G

function [ c ] = gen_matrix_c( m )
%GEN_COEFFS generate the coefficients involved in the B matrix
%   Detailed explanation goes here

c=zeros(m);
for i=1:m
    c(i,1)=1/(i+1);
end

for j=2:m
    for i=m:-1:2
        c(i-1,j)=c(i,j-1)*j/i;
    end
end

end