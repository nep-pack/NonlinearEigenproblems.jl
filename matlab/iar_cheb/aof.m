function [evps,VV,H]=aof(y0op,a,b,n,N,x0)
%  Arnoldi's method on a function (AOF)
%     [evps,VV,H]=aof(y0op,a,b,n,N,x0)
%     [evps,VV,H]=aof(y0op,a,b,n,N)
%   
%  Carries out N steps of method Arnoldi's method on a function
%  with Chebyshev scalar product on interval [a,b]. 
% 
%  Example: linear evp
%    n=80;    % construct large matrix with one large isolated eigenvalue
%    A0=randn(n); [Q,R]=qr(A0);
%    A0=Q'*diag([30,(randn(1,n-1)+1i*randn(1,n-1))])*Q;
%    cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
%    y0comp=@(X,Y) A0*(X*cheb_vect(0,X))-(Y*cheb_vect(0,Y))
%    [evps,V2,H,V]=aof(y0comp ,-1,1,n,10);
%    shouldbezero=min(min(evps)-1./eig(A0))  
%     
%  
%
%
    if (nargin==5)
        x0=ones(n,1);
    end



    orthtime=0; matvectime=0; sztime=0;

    x0=x0/norm(x0);
    H=[];
    V=zeros((N+1)*n,N+1);
    V(1:n,1)=x0;
    k0=1;
    for k=1:N
        % Apply matrix vector product to last column
        starttime=cputime;
        X=reshape(V(1:(n*k0),k),n,k0);
        Y=mat_vec_prod(X,n,y0op,a,b);
        y=reshape(Y,prod(size(Y)),1);
        matvectime=matvectime+(cputime-starttime);

        % orthogonalization and H-update
        starttime=cputime;
        [yo,h]=orthogonalize(V(1:((k0+1)*n),1:k),y,n,k,k0);
        orthtime=orthtime+(cputime-starttime);
        H(1:k+1,k)=h;
        starttime=cputime;
        V(1:length(yo),k+1)=yo; % save the orthognalized value
        kk0=[k,k0];
        % Only expand with new row if "numerically non-zero" 
        % This trick actually causes the convergence to stagnate
        % for some larger eigenvalues (k=40, sqrt-problem)
        %        if (norm(yo((end-n):end))>10*eps)
            k0=k0+1;
            %end
        sztime=sztime+(cputime-starttime);

    end
    
    % Print the profiling
    matvectime
    orthtime
    sztime
    

    % Compute the eigs of the Hessenberg matrix and Ritz vec's
    [vv,D]=eig(H(1:end-1,:));
    evps=1./diag(D);
    VV=V(1:n,1:end-1)*vv;
    VV=[VV;V(n+(1:n),1:end-1)*vv];
    
end 

function [w,h]=orthogonalize(V,w,n,kmax,k0)

    h=V'*w ;
    w=w-V*h ;
    g=V'*w ; % g is zero in exact arithmetic
    
    if (norm(g)>10*eps)  % only reorth. if necessary
        w=w-V*g ;
        h=h+g ;
    end
    
    h=[h;norm(w)] ;
    w=w/h(end) ; 
end

function Y=mat_vec_prod(X,n,y0op,a,b)
%  The abstract matrix vector product
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


