function [ Fval, Fder ] = FF( x, A )
%FF Summary of this function goes here
%   Detailed explanation goes here


Fval=x'*A*x;
Fder=x'*(A+A');

end

