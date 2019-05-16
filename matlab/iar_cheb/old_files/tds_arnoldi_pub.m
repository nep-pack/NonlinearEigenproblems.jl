function [evps,VV,H,V,iterdata]=tds_arnoldi_pub(tds,varargin)
%  tds_arnoldi: Compute eigenvalues of delay system with cheb. Arnoldi
%
%   [evps,VV]=tds_arnoldi_pub(tds)
%   [evps,VV]=tds_arnoldi_pub(tds,x0,N)
%  
%   Computes the eigenvalues of the time-delay system tds with the
%   delay Arnoldi scheme with starting vector x0 and N iterations.
%   
%   Example:  
%    A0=randn(3); A1=randn(3); A2=randn(3);
%    tds.A={A0,A1,A2};
%    tds.hA=[0,1,1.5];
%    [evps,VV]=tds_arnoldi_pub(tds);
%    s=min(evps);
%    M=-s*eye(3)+A0+A1*exp(-tds.hA(2)*s)+A2*exp(-tds.hA(3)*s);
%    shouldbezero=min(svd(M))
% 
%  Author: Elias Jarlebring 2010
%  Please cite:
%     * E. Jarlebring, M. Bennedich, G. Mele, E. Ringh, P. Upadhyaya, NEP-PACK: A Julia package for nonlinear eigenproblems, arxiv preprint, 2018
%  and
%     * E. Jarlebring, K. Meerbergen and W. Michiels, A Krylov method for the delay eigenvalue problem, SIAM J. Sci. Comput., 32(6):3278-3300, 2010
%
%
    


    n=length(tds.A{1});
    if (length(varargin)==0)
        x0=randn(n,1);
        N=20;
    else    
        x0=varargin{1}; 
        N=varargin{2};
    end

    if (size(x0,2)>1  | size(x0,1)~=length(tds.A{1}) | norm(x0)==0)
        error('x0 must be a nonzero vector of length n');
    end

    x0=x0/norm(x0); % normalize starting vector
    V=x0;
    H=[];

    orthtime=0; matvectime=0; % profiling

    taumax=max([tds.hA]);

    
    %   Construct the backward solve using lu-decomposition:
    %   Asuminvfun=@(x) Asum\x;

    Asum=tds.A{1}; 
    for k=2:length(tds.A)
        Asum=Asum+tds.A{k}; 
    end

    if (length(Asum)==1)        
        Asuminvfun=@(x) x/Asum;
    else
        if (issparse(Asum))
            [L,U,P,Q]=lu(Asum);
            Asuminvfun=@(x) Q*(U\(L\(P*x)));
            
            % For older version of octave:
            %[L,U,P]=lu(Asum);
            %Asuminvfun=@(x) U\(L\(P*x));
        else
            [L,U,P]=lu(Asum);
            Asuminvfun=@(x) U\(L\(P*x));
        end
    end

    %% Time for the Arnoldi iteration
    for k=1:N
        k0=size(V,1)/n; 
        
        % Apply matrix vector product to last column
        starttime=cputime;
        Y=reshape(V(:,end),n,k0);
        X=SigmaInvPiY(Y,tds,n,taumax,Asuminvfun);
        x=reshape(X,prod(size(X)),1);
        
        if (k==1)
           iterdata.X=X;
        end
        % orthogonalization and H-update
        V=[V;zeros(n,size(V,2))];         
        matvectime=matvectime+cputime-starttime;

        starttime=cputime;
        [xorth,h]=orthogonalize(V,x);
        H(1:k+1,k)=h;

        V=[V,xorth];            


        
        orthtime=orthtime+cputime-starttime;
        
    end
    
    % Profiling
    iterdata.matvectime=matvectime;
    iterdata.orthtime=orthtime; 
    
    % Only return eigenvalues which have a residual less than some
    % tolerance 
    [evps,VV]=fetch_good_ritzpairs(tds,n,H,V,1e-4); 

end 

function [w,h]=orthogonalize(v,w)
    h = v'*w ;
    w = w - v * h ;
    g = v'*w ;
    w = w - v * g ;
    h = h + g ;
    h = [h ; norm(w)] ;
    w = w / h(size(h,1)) ; 
end



function X=SigmaInvPiY(Yhat,tds,n,taumax,Asuminvfun)
%  The matrix vector product
%

    k=size(Yhat,2); 
    

    %  Construct the band matrix
    if (k==1)
        L=taumax*[0.5]; % Octave friendly 
    else
        dv=1:(k+1);
        L=diag(1./dv,0)+diag(-[1./dv(1:end-2)],2);
        L=L(1:end-1,:);
        L(1,1)=2;
        L=(taumax/4)*L;
        L=L(:,1:end-1);        
    end
    L=sparse(L);

    tic 
    %% Compute the "shifted down part"    
    Z=Yhat*L';


    %% Compute xhat
    ysum=sum(Yhat,2);
    % zsum=sum(Z,2);
    
    time1=toc;
    tic
    xx=ysum;
    for j=1:length(tds.hA)
        Tsum=zeros(size(xx));
        % Potential optimization: write as diagonal matrix product
        for i=1:size(Z,2)
            Tsum=Tsum+cheb_eval(1-2*tds.hA(j)/taumax,i)*Z(:,i);
        end
        xx=xx-tds.A{j}*Tsum;
    end

        
    xhat=Asuminvfun(xx);


    
    X=[xhat,Z];
    time2=toc;
    
end
    
function t=cheb_eval(x,n)    
%  Evaluate the chebyshev polynomial
    t=cos(n*acos(x));
end

function [ev,Vv]=fetch_good_ritzpairs(tds,n,H,V,TOL)
% Return the ritz pairs corresponding to approximations
% with a residual less than TOL.
    [vv,D]=eig(H(1:end-1,:));
    allevs=1./diag(D);
    
    % Try to extract some eigenvectors
    VV1=V(1:n,1:end-1)*vv;
    
    ev=[]; Vv=[];
    
    for k=1:size(VV1,2)
        x=VV1(:,k);
        x=x/norm(x);
        s=allevs(k);
        
        r=-s*x;
        for i=1:length(tds.A)
            r=r+tds.A{i}*(exp(-s*tds.hA(i))*x);
        end
        if (norm(r)<TOL)
            ev=[ev,s];
            Vv=[Vv,x];
        end
        
    end
end
