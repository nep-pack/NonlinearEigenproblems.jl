% The GUN example for the NLEVP toolbox
%
% The following script generates the figures for the example with a
% with a square square root terms. 
%
% 
% If you use this code, please cite the paper where the example was
% originally published:
%
% 
%
% 
% 


% Load problem data
load gun_orig_data;

% original NEP before scaling
M_org=@(l) A{1}-l*A{2}+1i*sqrt(l)*A{3}+1i*sqrt(l-108.8774^2)*A{4};

%% Carry out transformation with a shift and a scaling. 
if (1)  % favorable for bilanczos
    sigma=300^2;
    scale=(300^2-200^2)*10;
else    % favorable for iar    
    sigma=300^2;
    scale=(300^2-200^2);
end
f1=@(l) l;
f2=@(l) sqrt(l*scale+sigma);
f3=@(l) sqrt(l*scale+sigma-108.8774^2);
A{1}=A{1}-sigma*A{2};
A{2}=-scale*A{2};
A{3}=1i*A{3};
A{4}=1i*A{4};
% This is the transformed NEP
M=@(l) A{1}+A{2}*f1(l)+A{3}*f2(l)+A{4}*f3(l);


%%% Create a SMSF: Sum of products of matrices and scalar functions.
%%% Store derivatives in tcoeffs and matrices in Av:
NN=200;  % number of derivatives to compute a priori
Av=A;
% tcoeffs contains the derivatives of f1,f2,f3,f4
tcoeffs=zeros(length(Av),NN);
tcoeffs(1,1)=0;
tcoeffs(2,1)=1;
for k=1:NN
    %tcoeffs(3,k)=  scaled_sqrt_der(scale, sigma,k-1)/factorial(k-1);
   tcoeffs(3,k)=  scaled_sqrt_der(scale, sigma,k);
   %tcoeffs(4,k)=  scaled_sqrt_der(scale, sigma-108.8774^2,k-1)/factorial(k-1);
   tcoeffs(4,k)=  scaled_sqrt_der(scale, sigma-108.8774^2,k);
end
% Scale derivates with a factorial as a precomputation
tcoeffs_scaled=zeros(size(tcoeffs));
for k=1:NN
    tcoeffs_scaled(:,k)=tcoeffs(:,k)/factorial(k);
end
% Make the derivative matrix sparse in order to avoid 
% zero multiplication corresponding to linear and constant term
tcoeffs_scaled=sparse(tcoeffs_scaled); 

display('Doing LU');
[L,U,p,q,R]=lu(M(0),'vector');
q2(q) = 1:length(q);
LT = L';        UT = U';    

% Actions and operations involving M.
clear nep;
nep.M_solve     =@(x)      luS(x,L,U,p,q2,R);
nep.Mt_solve    =@(x)      R(p,:)\(LT\(UT\x(q)));
%nep.Mt_solve    =@(x)      M(0)'\x;
nep.M_lin_comb  =@(Y,i0)   smsf_lin_comb(Y,i0,Av,tcoeffs,0);
nep.Mt_lin_comb =@(Y,i0)   smsf_lin_comb(Y,i0,Av,tcoeffs,1);
nep.M           =          M;



randn('seed',0);
rand('seed',0);

v0 =randn(n,1);        
v0=v0/norm(v0);


u0=v0;
u0=u0/(u0'*nep.M_lin_comb(v0,1));


%% Run infinite bilanczos and IAR
% (and TIAR for m=20 for comparison) 

m=30; % number of iterations

tic
result_iar=iar(nep,v0,m);
iar_time=toc;

tic
nep.left_right_scalar_prod = @(nep,At,B,ma,mb) smsf_left_right_scalar_prod(nep,At,B,ma,mb,Av,tcoeffs,tcoeffs_scaled);
result_bilan=infbilanczos(nep,v0,u0,m);
bilanczos_time=toc;

tic
result_tiar=tiar(nep,v0,20);
tiar_time=toc;

H=result_bilan.T;
Hv={result_iar.H,H};
itertime_all={result_iar.timing_iteration,result_bilan.timing_iteration};
%toc
figure(1);clf;
figure(2);clf;
col='k--';
P=zeros(2,2);
for k=1:2
    H=Hv{k};
    for j=1:min(20,length(result_iar.eigvals))
        lstar=result_iar.eigvals(j);
        ev=[];
        for i=1:m-1
            eig0=eig(full(H(1:i,1:i)));
            ll=1./eig0;
            e=min(abs(lstar-ll));
            ev=[ev,e];
        end
        figure(1);
        P(k,1)=semilogy(itertime_all{k}(1:length(ev)),ev,col);
        hold on;
        xlabel('wall time')
        ylabel('eigenvalue error')
        figure(2);
        P(k,2)=semilogy(ev,col)        
        hold on;
        xlabel('iteration')
        ylabel('eigenvalue error')

    end
    col='r-*';
end
l=legend([P(2,1), P(1,1)],...
         ['Infinite bi-Lanczos'],...
         ['IAR']);

l=legend([P(2,2), P(1,2)],...
         ['Infinite bi-Lanczos'],...
         ['IAR']);


fprintf('Total time for IAR: %f\n',result_iar.timing_iteration(end))
fprintf('Total time for infinite bilanczos: %f\n',result_bilan.timing_iteration(end))
fprintf('Scalar product time: %f\n',sum(result_bilan.timing_scalarprod))