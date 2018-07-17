function [ D ] = BB_precomputation( nep, m_max )
%B_precomputation Precomputation in Infinite Lanczos
%   The matrix B can be block diagonalized. This function precomputes the 
%   block diagonal elements. 
%
%   The method is based on the idea that B can be block diagonalized as
%   explained in
%   Kaveh, Ali, and Hossein Rahami:  
%   "Block circulant matrices and applications in free vibration analysis 
%   of cyclically repetitive structures." 
%   Acta Mechanica 217.1 (2011): 51-62.

% precompute block-diagonals
[alpha,beta] = G(m_max);
D=cell(2*m_max);
for j=1:2*m_max
    D{j}=@(v,m) alpha(j,m)*v+beta(j,m)*(nep.A1*v);
end

end