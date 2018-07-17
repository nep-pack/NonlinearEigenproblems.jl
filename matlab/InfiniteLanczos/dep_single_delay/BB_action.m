function [ Z ] = BB_action( nep, Y )
%BB_ACTION Compute the product B*Y

Z=0*Y;
for j=1:size(Y,2)
    Z(:,j)=BB_action_vec( nep, Y(:,j) );
end


end


function [ z ] = BB_action_vec( nep, y )
%BB_ACTION Compute the product B*y 
%   Exploit the Block-Circulant structure of B in order to efficiently
%   compute the matrix-vector product B*y
%
%   The method is based on the idea that B can be block diagonalized as
%   explained in
%   Kaveh, Ali, and Hossein Rahami:  
%   "Block circulant matrices and applications in free vibration analysis 
%   of cyclically repetitive structures." 
%   Acta Mechanica 217.1 (2011): 51-62.


n=nep.n;                % size of the problem
m=length(y)/n;          % number of blocks
%pause
alpha=nep.alpha; beta=nep.beta;

y=reshape(y,n,m);
y=fliplr(y);
y=[y zeros(n,m)];
y=reshape(y,n,2*m);
z=(ifft(y')*sqrt(2*m))';
%fprintf("FFT on a matrix of size (%d,%d)\n",size(y',1),size(y',2));


%A1z=nep.A1*z;
%for j=1:2*m
%    z(:,j)=alpha(j,m)*z(:,j)+beta(j,m)*A1z(:,j);
%end
%z=z*alpha(1:2*m)+(nep.A1*z)*beta(1:2*m,:);
z=bsxfun(@times,z,alpha(1:2*m,m).')+bsxfun(@times,nep.A1*z,beta(1:2*m,m).');


z=(fft(z')/sqrt(2*m))'; 
z=z(:);   
z=z(1:m*n);


end
