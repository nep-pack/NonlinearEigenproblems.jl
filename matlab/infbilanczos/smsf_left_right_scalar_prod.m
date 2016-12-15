function c=smsf_left_right_scalar_prod(nep,At,B,ma,mb,Av,tcoeffs,tcoeffs_scaled)
%  computes a left-right scalar product specifically for the
%  structure "sum of matrices scalar functions" 
%
%  In order to use it set:
%  nep.left_right_scalar_prod=@(nep,At,B,ma,mb) sqrt_left_right_scalar_prod(nep,At,B,ma,mb,Av,tcoeffs,tcoeffs_scaled);

persistent Av_projected;
if (~exist('Av_projected'))  % trick to allocate it once 
    Av_projected=cell(size(Av)); 
    mm=200;
    for i=1:length(Av) 
        Av_projected{i}=zeros(mm,mm);
    end
end
m0=min(ma,mb);
for i=1:length(Av)  
    Av_projected{i}(1:ma,1:mb)=At(:,1:m0)'*(Av{i}*B(:,1:mb));
end
c=0;
for k=1:length(Av)
    II=reshape(ones(ma,1)*(1:mb)+(1:ma)'*ones(1,mb),ma*mb,1);
    T=tcoeffs_scaled(k,II-1);
    c=c+T*reshape(Av_projected{k}(1:ma,1:mb),ma*mb,1);
end

%% Didactic version of above  (but slower)
%for k=1:length(Av)
%    for l=1:mb
%        for j=1:ma
%           c=c+Av_projected{k}(j,l)*(tcoeffs(k,j+l-1)/factorial(j+l-1));
%        end
%    end
%end

c=-c;


