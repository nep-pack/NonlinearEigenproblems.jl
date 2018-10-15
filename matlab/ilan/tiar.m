function [ W, E ] = tiar( nep, k, opts )
%tiar Tensor Infinite Arnoldi for nonlinear eigenvalue problems
%
%   [W,E] = tiar(nep) produces a column vector E of nep's 6 smallest
%   magnitude eigenvalues and a full matrix W whose columns are the 
%   corresponding eigenvectors so that M(E(i))*W(:,i) = 0. nep is a 
%   structure describing the nonlinear eigenvalue problem.
%
%   nep.n:   size of problem
%
%   nep.M:   function handle for the evaluation of M(l)*v where
%            (l[scalar], v[vector])
%                           [ @(l,v) {vector M(l)*v} ]
%
%   (optional*): nep.Md:  function handle for the evaluation of M^(j)(0)v 
%                where M^(j)(0) is the j-th derivative of M(l) 
%                and (j[scalar], v[vector]). nep.Md or nep.Mlincomb has
%                to be provided.
%                           [ @(j,v) {vector M^(j)(0)*v} ]
%
%   (optional*): nep.Mlincomb: function handle for the evaluation of
%                M^(1)(0)X(:,1)+M^(2)(0)X(:,2)+...+M^(j)(0) X(:,j)
%                where M^(s)(0) is the s-th derivative of M(l) 
%                and (j[scalar], v[vector]). nep.Md or nep.Mlincomb has
%                to be provided.
%                           [ @(X) {vector \sum_{s=1}^j M^(s)(0)X(:,s)} ]
%
%   nep.Minv: function handle for the solution of linear systems 
%                 M(0)x=v where (v[vector])
%                            [ @(v) {solution of M(0)x=v } ]
%
%   (optional) nep.resnorm: function handle for residual norm where
%                           (l[scalar], v[vector])
%                            [ @(l,v) {estimated residual norm} ]
%
%   tiar(nep,k) produces a column vector E of nep's K smallest magnitude
%   eigenvalues and a full matrix W whose columns are the corresponding 
%   eigenvectors so that M(E(i))*W(:,i) = 0
%
%   tiar(A,k,opts) specify options:
%   opts.maxit: maximum number of iterations [integer | {100}]
%   opts.tol: convergence: resnorm < tol [scalar | {100*eps}]
%   opts.v0: starting vector [n-by-1 vector | {randomly generated}]
%   opts.sigma: compute the eigenvalues close to the shift sigma. If the 
%               shift is specified, nep.Md (if provided) has to be the 
%               function handle for the evaluation of M^(j)(sigma)v,  
%               nep.Mlincomb (if provided) has to be the 
%               function handle for the evaluation of 
%               M^(1)(sigma)X(:,1)+M^(2)(sigma)X(:,2)+...+M^(j)(sigma) X(:,j)
%               and nep.Minv the function handle for the solution 
%               of linear systems M(sigma)x=v [scalar |{0}]
%   opts.gamma: rescale the problem
%   opts.p: the convergence is checked every p iterations [integer |{1}]
%   opts.disp: diagnostic information display level [{0} | 1 ]   
%
%   Example (Quadratic Eigenvalue Problem):
%   compute eigenpairs of M(l)=l^2*A2+l*A1+A0
%
%   clear nep
%   n=1000;
%   A0=sprand(n,n,1/n)+speye(n);
%   A1=sprand(n,n,1/n);
%   A2=sprand(n,n,1/n);
%   sigma=1+1i; gamma=3;
%   nep.M =@(l,v) l^2*A2*v+l*A1*v+A0*v;
%   nep.Md=@(j,v) (j==1)*(A1*v+2*sigma*A2*v)+(j==2)*2*A2*v;
%   [L,U,P,Q] = lu(A0+sigma*A1+sigma^2*A2);
%   nep.Minv=@(v) Q*(U\(L\(P*v)));
%   nep.n=n;
%   nep.resnorm=@(l,v) norm(nep.M(l,v))/((abs(l)^2*norm(A2,inf)+abs(l)*norm(A1,inf)+norm(A0,inf))*norm(v));
%   opts.maxit=100; opts.tol=1e-10; opts.disp=1;  k=20;
%   opts.sigma=sigma;    opts.gamma=gamma;  opts.p=7;
%   [ W, E ] = tiar( nep, k, opts );
%   figure; plot(real(E),imag(E),'k*');
%   hold on; plot(real(sigma),imag(sigma),'ro')
%   xlabel('real'); ylabel('imag')
%   title('Converged Eigenvalues')
%
%   The following example shows that, for some problems, the algorithm is
%   more efficient if nep.Mlincomb is provided. Test it 
%
%   Example (Delay Eigenvalue Problem):
%   compute eigenpairs of M(l)=-l*I+A0+exp(-l)*A1
%
%   clear nep
%   n=1000;
%   A0 = spdiags(rand(n,3), -1:1, n, n);
%   A1 = spdiags(rand(n,3), -1:1, n, n);
%   nep.M =@(l,v)  -l*v+A0*v+exp(-l)*(A1*v);
% 
%   AA=cell(3,1);
%   AA{1}=@(v) v+A1*v;
%   AA{2}=@(v) A1*v;
%   AA{3}=@(v) A1*v;
%   nep.Md=@(j,v) AA{min(j,3)}(v)*(-1)^j;
%   [L,U,P,Q] = lu(A0+A1);
% 
%   nep.Minv=@(v) Q*(U\(L\(P*v)));
%   nep.n=n;
%   nep.resnorm=@(l,v) norm(nep.M(l,v))/((abs(l)+norm(A0,inf)+abs(exp(-l))*norm(A1,inf))*norm(v));
%   opts.maxit=50; opts.tol=1e-10; opts.disp=0; opts.v0=ones(n,1); opts.p=inf; k=inf;
% 
%   tic; 
%   [ W1, E1 ] = tiar( nep, k, opts ); 
%   toc;
%  
%   tic; 
%   nep.Mlincomb=@(X,j) -X(:,1)+A1*(sum(bsxfun(@times, X(:,1:j),(-1).^(1:j)),2));
%   [ W2, E2 ] = tiar( nep, k, opts );
%   toc;
%
%   If you use this code, please cite the following article:
%   E. Jarlebring, G. Mele, and O. Runborg
%   The waveguide eigenvalue problem and the tensor infinite Arnoldi method
%   SIAM J. Sci. Comput., accepted for publication, 2017
%
%   The code is published to improve reproducability, and is free to use.
%   However, the code comes with no waranty an no guarantee. The authors 
%   has no responsability for outcomes from the ussage of this code
%   Use at your own risk.
%
%   author:  Elias Jarlebring and Giampaolo Mele
%   version: 25/05/2017

n        =  nep.n;
Minv     =  nep.Minv;
M        =  nep.M;

% check input
if isfield(nep,'resnorm')
    resnorm = nep.resnorm;
else
    resnorm = @(l,v) norm(M(l,v));
end 

if ~exist('k','var')
    k=6;
end
if (isfield(nep,'Mlincomb'))
    Mlincomb=nep.Mlincomb;
elseif isfield(nep,'Md')
    % Compute linear combination of derivatives
    % in the naive way by calling Md several times:
    Mlincomb=@(X,j) sum(bsxfun(@(jj,y) nep.Md(jj,X(:,jj)), 1:j, zeros(n,1)),2);
else    
    error('You must provide either Md or Mlincomb');
end



if ~exist('opts','var')
    maxit=100;  tol=100*eps;
    disp=0;     v0=rand(n,1);
    sigma=0;    gamma=1;
else
    
    if isfield(opts,'maxit')
        maxit        =  opts.maxit;
    else 
        maxit=100;
    end
    
    
    if isfield(opts,'tol')
        tol      =  opts.tol;
    else 
        tol=1e-12;
    end
    
    if isfield(opts,'disp')
        disp      =  opts.disp;
    else 
        disp=0;
    end
    
    if isfield(opts,'v0')
        v0      =  opts.v0;
    else 
        v0=rand(n,1);
    end
    
    if isfield(opts,'sigma')
        sigma=opts.sigma;
    else 
        sigma=0;
    end
    
    if isfield(opts,'gamma')
        gamma=opts.gamma;
    else 
        gamma=1;
    end
    
    if isfield(opts,'p')
        p=opts.p;
        if p>maxit
            p=maxit;
        end
    else 
        p=1;
    end
    
end
    
% initialize tiar factorization
Z=zeros(n,maxit+1);     H=zeros(maxit+1,maxit);
size(v0)
Z(:,1)=v0;              Z(:,1)=Z(:,1)/norm(Z(:,1));
a(1,1,1)=1;             W=[];   E=[];

% expand the tiar factorization
j=1; conv_eig=0;    kp=1;
while (j<=maxit-1)&&(conv_eig<k)
  
    if disp==1
        if p==1
            fprintf('Current iteration: %d, converged eigenpairs: %d  \r',j,min(conv_eig,k));
        else
            fprintf('Current iteration: %d, converged eigenpairs (at iteration %d): %d  \r',j,kp,min(conv_eig,k));
        end
    end
    
    % computation of y(:,2), ..., y(:,k+1)     
    AK(1:j,1:j)=a(1:j,j,1:j); 
    y(:,2:j+1) = bsxfun (@rdivide, Z(:,1:j)*AK.', 1:j);

    % computation of y(:,1)
    if (gamma~=1)
        ZZ=bsxfun(@times, y(:,1:(j+1)), gamma.^(0:j));    
    else
        ZZ=y(:,1:(j+1));
    end
    ZZ(:,1)=0;
    y(:,1)=Mlincomb(ones(j,1),ZZ(:,2:end));
    y(:,1)=-Minv(y(:,1));
    
    % Gram–Schmidt orthogonalization in Z
    t=zeros(j,1);    t(1:j)=Z(:,1:j)'*y(:,1);
    Z(:,j+1)=y(:,1)-Z(:,1:j)*t(1:j);   
    
    % Gram–Schmidt re-orthogonalization in Z
    tt=zeros(j,1);  tt(1:j)=Z(:,1:j)'*Z(:,j+1);
    Z(:,j+1)=Z(:,j+1)-Z(:,1:j)*tt(1:j);   
    t=t+tt;         t(j+1)=norm(Z(:,j+1));
    Z(:,j+1)=Z(:,j+1)/t(j+1);
    
    % extension of the tiar factorization 
    a(j+1,:,1:j+1)=0;   
    
    % compute the matrix G
    for l=1:j+1
        for i=2:j+1
            g(i,l)=a(i-1,j,l)/(i-1);
        end
        g(1,l)=t(l);    g(i,j+1)=0;
    end
    
    % compute h (orthogonalization with tensors factorization)
    h=zeros(j,1);
    for l=1:j
        Al=a(1:j,1:j,l);
        h=h+Al'*g(1:j,l);
    end
    
    % compute the matrix F
    for l=1:j
        Al=a(1:j+1,1:j,l);
        f(1:j+1,l)=g(1:j+1,l)-Al*h;
    end
    
    for i=1:j+1
        f(i,j+1)=g(i,j+1);
    end
            
    % re-orthogonalization
    % compute hh (re-orthogonalization with tensors factorization)
    hh=zeros(j,1);
    for l=1:j
        Al=a(1:j,1:j,l);
        hh=hh+Al'*f(1:j,l);
    end
    
    % compute the matrix FF
    for l=1:j
        Al=a(1:j+1,1:j,l);
        ff(1:j+1,l)=f(1:j+1,l)-Al*hh;
    end
    
    for i=1:j+1
        ff(i,j+1)=f(i,j+1);
    end
    
    % update the orthogonalization coefficients
    h=h+hh; f=ff;
    
    beta=norm(f,'fro');
    
    % extend the matrix H
    H(1:j,j)=h; H(j+1,j)=beta;
    
    % extend the tensor
    for i=1:j+1
        for l=1:j+1
            a(i,j+1,l)=f(i,l)/beta;
        end
    end

    % compute Ritz pairs every P iteration
    if (rem(j,p)==0)||(j==maxit-1)
        kp=k;

        % partial extraction of the basis
        AA(1:j+1,1:j+1,1)=a(1,:,:);  AA=AA.';
        V=Z(:,1:j+1)*AA;

        [RitzVec, RitzVal] = eig(H(1:j,1:j));
        RitzVal=diag(RitzVal);      RitzVal=sigma+gamma./RitzVal;
        RitzVec=V(:,1:j)*RitzVec;   RitzVec=RitzVec(1:n,:);

        conv_eig=0;
        for i=1:j
            RitzVec(:,i)=RitzVec(:,i)/norm(RitzVec(:,i));
            ERR(i,j)=resnorm(RitzVal(i),RitzVec(:,i));
            % count converged Ritz pairs
            if ERR(i,j)<tol
                conv_eig=conv_eig+1;
            end
        end
        % sort the error vector
        [ERR(1:j,j),I]=sort(ERR(1:j,j));

        % return the converged Ritz pairs
        if(conv_eig>=k)||(j==maxit-1)   
            W=RitzVec(:,I(1:min(k,conv_eig))); 
            E=RitzVal(I(1:min(k,conv_eig)));
        end
        
    end
    
    j=j+1;
end

if disp==1
    fprintf('\rTotal number of iterations: %d\nConverged eigenpairs: %d \n',j,min(conv_eig,k));
    if (p<=j/2)&&(p<j-1)
        figure
        for i=1:size(ERR,2)   
            semilogy(p:p:j-1,ERR(i,p:p:j-1),'-k');   
            hold on
        end
        axis([p j eps/100 1])
        xlabel('Iteration'); ylabel('resnorm'); title('Convergence history');
    end
end

if conv_eig<k
    warning off backtrace;
    warning('\rOnly %d of %d eigenvalues computed.\nIncrease opts.maxit and/or decrease opts.tol',conv_eig,k);
end



end
