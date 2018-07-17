function [ c ] = gen_coeffs( m )
%GEN_COEFFS generate the coefficients involved in the B matrix
%   Detailed explanation goes here


c=zeros(m);
for i=1:m
    c(i,1)=1/i;
end

for j=1:m-1
    for i=m:-1:2
        c(i-1,j+1)=c(i,j)*j/(i-1);
    end
end

end

