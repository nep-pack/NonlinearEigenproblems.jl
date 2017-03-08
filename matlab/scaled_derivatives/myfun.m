function [ f ] = myfun( x, k )
%MYFUN auxiliary function made to evaluate funm
%   In order to use funm you need to compute the derivatives

% Return kth derivative of exp + cos at X.
g = mod(ceil(k/2),2);
if mod(k,2)
    f = exp(x) + sin(x)*(-1)^g;
else
    f = exp(x) + cos(x)*(-1)^g;
end


end

