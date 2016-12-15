function z=smsf_lin_comb(Y,i0,Av,tcoeffs,do_transpose)
% computes the linear combination corresponding to Av and tcoeffs 
%

n=length(Av{1});
z=zeros(n,1);
if (i0<1)
    error('Not implemented for i0<1'); 
end
for i=1:length(Av)
%% Equivalent but slower
%    for j=1:size(Y,2)
%        if (do_transpose)
%            z=z+Av{i}'*(Y(:,j)*tcoeffs(i,i0+j-1)');
%        else
%            z=z+Av{i}*(Y(:,j)*tcoeffs(i,i0+j-1));         
%        end
%    end
     if (do_transpose)
         z=z+Av{i}'*(Y*tcoeffs(i,i0:(i0+size(Y,2)-1))');
     else
         z=z+Av{i}*(Y*tcoeffs(i,i0:(i0+size(Y,2)-1)).');         
     end
end
