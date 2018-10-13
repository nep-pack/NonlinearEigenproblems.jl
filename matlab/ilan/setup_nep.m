function nep=setup_nep(f,A,mu,m,fD)
    % setup_nep generate all the needed data for solving the nep
    % f are the nonlinearities
    % A are the matrices
    % mu is the shift
    % m is the number of derivatives wanted
    %
    % This function will create an interface

    % THE PROBLEM HAS TO BE IN SPMF FORMAT
    % M(l) = f_1(l)A_1+...+f_p(l)A_p

    p=length(f);
    
    % precompute all the derivatives
    % precompute all the derivatives
    if isempty(fD)
        fD=zeros(m,p);
        for t=1:p
            fD(:,t)=compute_fder(f{t},mu,m);
        end
    end
        
    % precompute the sparse lu factorization for
    % M0
    M0=compute_Mder_eval(0,A,fD);
    if issparse(M0)
        [L,U,P,Q,R] = lu(M0);
        Minv=@(v) Q*(U\(L\(P*(R\v))));
    else
        [L,U,P] = lu(M0);
        Minv=@(v) U\(L\(P*v));
    end

    % precompute the norms of the matrix
    % coefficients
    nA=zeros(p);
    for t=1:p
        nA(t)=norm(A{t},'inf');
    end
    
    % precompute Hankel derivative matrices
    FHD=cell(length(f),1);
    for t=1:length(f)
        FHD{t}=Fmat(fD,t,ceil((m-2)/2));
    end
    
    
    % store in a structure
    nep.n=size(A{1},1); nep.f=f; nep.A=A; nep.fD=fD;
    nep.Mlincomb=@(alpha,V) compute_Mlincomb(alpha,V,A,fD);
    nep.Mlincombx=@(alpha,x) compute_Mlincomb_x(alpha,x,A,fD);
    nep.Mdeval=@(j) compute_Mder_eval(j,A,fD);
    nep.Md=@(j,v) compute_Mder(j,A,fD,v);
    nep.Minv=@(v) Minv(v);
    nep.resnorm=@(l,v) compute_resnorm(A,nA,f,l,v);
    nep.Meval=@(l) compute_Meval(A,f,l);
    nep.M=@(l,v) compute_M(A,f,l,v);
    nep.FHD=FHD;
end

function z=compute_Mlincomb(alpha,V,A,fD)
    % compute
    % z=alpha(1)M^(1)v1+...+alpha(k)M^(k)vk
    n=size(V,1); p=size(fD,2); k=length(alpha);
    z=zeros(n,1);
    fDa=bsxfun(@times,fD(2:k+1,:),alpha);
    for t=1:p
        z=z+A{t}*(V*fDa(:,t));
    end
end

function z=compute_Mlincomb_x(alpha,x,A,fD)
    % compute
    % z=alpha(1)M^(1)x+...+alpha(k)M^(k)x
    n=size(x,1); p=size(fD,2);
    k=length(alpha);
    z=zeros(n,1);
    fDa=bsxfun(@times,fD(2:k+1,:),alpha);
    for t=1:p
        z=z+sum(fDa(1:k,t))*(A{t}*x);
    end
end

function Mj=compute_Mder_eval(j,A,fD)
    % compute M^(j)
    p=size(fD,2);
    n=size(A{1},1);
    Mj=sparse(n,n);
    for t=1:p
        Mj=Mj+fD(j+1,t)*A{t};
    end
end

function z=compute_Mder(j,A,fD,v)
    % compute M^(j)*v
    p=size(fD,2);
    n=size(A{1},1);
    z=zeros(n,1);
    for t=1:p
        z=z+fD(j+1,t)*(A{t}*v);
    end
end

function fD = compute_fder(f,mu,m)
%COMPUTE_FDER compute derivatives of a function in mu
%   Compute the first m detivatives, evaluated in mu, of the function f
%   fD(j)=f^(j)(mu)
    S=mu*diag(ones(m,1))+diag(1:m-1,-1); 
    fD=f(S); 
    fD=fD(:,1);
end

function err = compute_resnorm(A,nA,f,l,v)
% Use the classic formula for computing the
% relative error for an eigenpair (l,v)
    p=length(f);
    n=size(A{1},1);
    v=v/norm(v);
    ns=0;
    x=zeros(n,1);
    for t=1:p
        x=x+f{t}(l)*(A{t}*v);       % sum over the nep nonlinearities
        ns=ns+nA(t)*abs(f{t}(l));   % scale factor
    end
    err=norm(x)/ns;
end

function M = compute_Meval(A,f,l)
% Use the classic formula for computing the
% relative error for an eigenpair (l,v)
    p=length(f);
    n=size(A{1},1);
    M=sparse(n,n);
    for t=1:p
        M=M+f{t}(l)*A{t};       % sum over the nep nonlinearities
    end
end

function z = compute_M(A,f,l,v)
% Use the classic formula for computing the
% relative error for an eigenpair (l,v)
    p=length(f);
    n=size(A{1},1);
    z=zeros(n,1);
    for t=1:p
        z=z+f{t}(l)*(A{t}*v);       % sum over the nep nonlinearities
    end
end


function F = Fmat(fD,t,k)
%FMAT Contruct the derivative-coefficients Hankel matrix
%   Construct the matrix
%  \begin{pmatrix}
%  f_t^{(1)} (\mu)      & f_t^{(2)}(\mu)	& f_t^{(3)}(\mu)	& \dots 	& f_t^{(k+1)}(\mu)	\\
%  f_t^{(2)} (\mu)      & f_t^{(3)}(\mu)	& \dots			& 	 	& f_t^{(k+2)}(\mu)	\\
%  f^{(3)} (\mu)		&  \dots		&			&		& f_t^{(k+3)}(\mu)	\\
%   \vdots              &			&			&		& \vdots		\\
%  f_t^{(k+1)} (\mu)	&  f_t^{(k+2)}(\mu)	& f_t^{(k+3)}(\mu)	& \dots		& f_t^{(2k+1)}(\mu)
%  \end{pmatrix}

F=zeros(k,k);
for i=1:k
    for j=1:k
        F(i,j)=fD(i+j,t);
    end
end


end