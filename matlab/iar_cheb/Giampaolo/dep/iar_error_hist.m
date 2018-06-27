function [ err, conv_eig, conv_vec ] = iar_error_hist( nep, V, H, mplot )
%RK_ERROR_HIST Summary of this function goes here
%   Detailed explanation goes here

conv_eig=[];
conv_vec=[];
m=size(V,2);
err=NaN(m,m);

for k=1:m-1    
    %k
%    [Z,S]=eig(H(1:k,1:k));
    [Z,S]=eig(H(1:k,1:k));
    
    for i=1:k
        s=1/S(i,i);
        z=V(:,1:k)*Z(:,i);  
        err(i,k)=nep.err(s,z(1:nep.n));
%        err(i,k)=min(svd(nep.MMeval(s)));
        
        % save eigenvalues if converged
        if err(i,k)<1e-5
            if isempty(conv_eig) || min(abs(conv_eig-s))>1e-6
                conv_eig=[conv_eig s];
                conv_vec=[conv_vec z(1:nep.n)];
            end                
        end
        
        
    end
    err(:,k)=sort(err(:,k),'ascend');
end


for i=1:m
    semilogy(err(i,:),mplot)
    hold on
end
ylim([1e-16 1e5])

end

