function [evps,V2,H,V]=aof_abs_cheb(y0op,a,b,n,N)
%   arnoldi's method on a function
%   y0op(X,Y)  
%   
%  
%   Example: Linear eigenvalue problem
%    n=80;    % (large matrix with one large isolated eigenvalue
%    A0=randn(n); [Q,R]=qr(A0);
%    A0=Q'*diag([30,(randn(1,n-1)+1i*randn(1,n-1))])*Q;
%    cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
%    y0comp=@(X,Y) A0*(X*cheb_vect(0,X))-(Y*cheb_vect(0,Y))
%    [evps,V2,H,V]=aof_abs(y0comp ,-1,1,n,10);
%    shouldbezero=min(min(evps)-1./eig(A0))  
%
%   Example: QEP
%     n=30;
%     M=randn(n); C=randn(n); K=randn(n);
%     a=-8; b=1; k=2/(b-a) ; c=(a+b)/(a-b);
%     cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
%     cheb2_vect=@(t,Z) sin( (1:(size(Z,2))) *acos(t))'./sin(acos(t));
%     Mterm=@(t,X) X(:,2:end)*(k*(1:(size(X,2)-1))'.*cheb2_vect(t,X(:,1:end-1)))
%     y0comp=@(X,Y) -K\(M*Mterm(c,X)+C*X*cheb_vect(c,X)+K*Y*cheb_vect(c,Y));
%     [evps,V2,H,V]=aof_abs(y0comp ,a,b,n,20);
%     shouldbezero=min(polyeig(K,C,M))-min(evps)
%     l=min(evps);
%     shouldbezero=min(eig(M*l^2+C*l+K))
%
    
    
    x0=ones(n,1);
    
    %    mat_vec_prod(randn(n,3),n,y0op,a,b)
    %    [Y0]=mat_vec_prod(x0,n,y0op,a,b)    
    %    keyboard
    [H,VV]=abs_arnoldi(@(X) mat_vec_prod(X,n,y0op,a,b),...
                      x0,N,@myscalarprod,@myaddscaleop);

    % Compute the Ritz values 
    [V,D]=eig(H(1:N,1:N));
    evps=1./diag(D);

    % Construct the Ritzvectors in V2
    V2=cell(length(V),1);
    for k=1:length(evps)
        v=V(:,k);
        V2{k}=0;
        for j=1:length(v)
            VVjfull=full(VV{j});
            V2{k}=myaddscaleop(v(j),VVjfull,1,V2{k});
        end
    end
    
end 


function Y=mat_vec_prod(X,n,y0op,a,b)
%%  The abstract matrix vector product
%

    k=size(X,2); 

    % Construct the band matrix:
    if (k==1)
        L=(b-a)/2;
    else
        L=diag(ones(k,1))+diag(-ones(k-2,1),2);
        L(1,1)=2;  % "Fix" assymetry
        L=((b-a)/4)*(diag(1:k)\L);
    end

    %% Compute the "shifted down part"    
    Y=X*L';

    %% Compute new element
    y0=y0op(X,[zeros(n,1),Y]);
    
    Y=[y0,Y];
end
    
function Z=myaddscaleop(c1,X1,c2,X2)
%% The abstract way of adding and scaling vectors. 
    n=size(X1,1);
    
    if ((norm(X1,1)==0) || (c1==0))
        Z=c2*X2; return;
    end

    if ((norm(X2,1)==0)  || (c2==0))
        Z=c1*X1; return;
    end

    % Dynamically increase the size
    if (size(X1,2)>=size(X2,2))
        Z=c1*X1;
        Z(:,1:size(X2,2))=Z(:,1:size(X2,2))+c2*X2;
    else
        Z=c2*X2;
        Z(:,1:size(X1,2))=Z(:,1:size(X1,2))+c1*X1;        
    end
end

% 
function v=myscalarprod(X,Y)
    n=size(X,1);
    if (size(X,1) ==0 || size(Y,1)==0)
        v=0; return;
    end

    i=min([size(X,2),size(Y,2)]);
    v=(reshape(X(:,1:i),i*n,1))'*(reshape(Y(:,1:i),i*n,1));
    %    v=sum(diag(X(:,1:i)'*Y(:,1:i)));
    
end


function x=eval_cheb_vect(X,t)
% Evaluate a function represented as a Chebyshev expansion 
% with coefficients given by the columns of X.
    x=zeros(size(X,1),1);
    for k=1:size(X,2)
        x=x+X(:,k)*cheb_eval(t,k-1);
    end
end