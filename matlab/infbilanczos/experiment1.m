% Experiment 1: A large QEP-DEP 
%
% The following script generates the figures for the example with a
% quadratic delay eigenvalue problem with random coefficient matrices
% 

clear all; close all;
% Problem data
randn('seed',1); rand('seed',0);        % need to reset both for sprandn
n=1000;
A0=sprandn(n,n,0.01);
A1=sprandn(n,n,0.01);                   
tau=1;
M=@(l) -l^2*speye(n)+A0+exp(-tau*l)*A1;
Mp=@(l) -2*l*speye(n)-tau*exp(-tau*l)*A1; % derivative of M. Needed
                                          % for Newton

%%% Create a SMSF: Sum of products of matrices and scalar functions.
%%% Store derivatives in tcoeffs and matrices in Av:
NN=200; % number of derivatives to compute a priori
Av={speye(n),A0,A1};
tcoeffs=zeros(3,NN);
tcoeffs(1,1)=0;
tcoeffs(1,2)=-2;
for i=1:NN
    tcoeffs(3,i)=(-tau)^i;
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
% For the operation M(0)^{-1} and M(0)^{-*} we use an LU.
[L,U,p,q,R]=lu(A0+A1,'vector');
q2(q) = 1:length(q);
LT = L';        UT = U';    
A0T = A0';      A1T = A1';

% Define the nonlinear eigenvalue problem NEP: Actions and operations involving M.
clear nep;

nep.M_solve     =@(x)           luS(x,L,U,p,q2,R);
nep.Mt_solve    =@(x)           R(p,:)\(LT\(UT\x(q)));
nep.M_lin_comb  =@(Y,i0)        smsf_lin_comb(Y,i0,Av,tcoeffs,0);
% nep.M_lin_comb  =@(Y,i0)        qep_dep_lin_comb(Y,A0,A1,tau,i0);
nep.Mt_lin_comb =@(Y,i0)        smsf_lin_comb(Y,i0,Av,tcoeffs,1);  
% nep.Mt_lin_comb =@(Y,i0)        qep_dep_lin_comb(Y,A0',A1',tau,i0);
nep.M           =               M;
nep.Mp          =               Mp;


%% Use specialized left-right scalar product
% 
% By commenting the following line out, you active the "generic" left-right-scalar-product.
nep.left_right_scalar_prod=@(nep,At,B,ma,mb) smsf_left_right_scalar_prod(nep,At,B,ma,mb,Av,tcoeffs,tcoeffs_scaled);
% nep.left_right_scalar_prod=@(nep,At,B,ma,mb) qep_dep_left_right_scalar_prod(nep,At,B,ma,mb,A0,A1,tau)

%% Functions to measure the backward error and codition number of the eigenvalues.
normalize =@(v) v/norm(v,1);
v0 =randn(n,1);         

normA0 = norm(full(A0),2);
normA1 = norm(full(A1),2);
nep.resnorm_l   = @(l)      norm(M(l)*normalize(M(l)\v0),1)/(abs(l)^2+norm(A0,1)+norm(A1,1)*abs(exp(-tau*l)));
nep.resnorm_lx  = @(l,x)    norm(M(l)*x,1) / ( abs(l)^2+norm(A0,1)+norm(A1,1)*abs(exp(-tau*l)) * norm(x,1));
nep.resnorm_ly  = @(l,y)    norm(y'*M(l),1) / (abs(l)^2+norm(A0,1)+norm(A1,1)*abs(exp(-tau*l)) * norm(y,1));
nep.resnorm_lxy = @(l,x,y)  max( norm(M(l)*x,1) , norm(y'*M(l),1) ) / ...
                                ( norm(x,1)*norm(y,1)*(abs(l)^2+norm(A0,1)+norm(A1,1)*abs(exp(-tau*l))) ); 
nep.cond_lxy    = @(l,x,y)  (abs(l^2)+normA0+abs(exp(-tau*l))*normA1) * norm(x,2) * norm(y,2) / ...
                                    (abs(l)*abs(y'*((-2*l*speye(n)-tau*exp(-tau*l)*A1)*x)));
nep.cond2_lxy    = @(l,x,y) ((abs(l^2)+normest(A0)+abs(exp(-tau*l))*normest(A1))) * norm(x,2) * norm(y,2) / ...
                                      (abs(l)*abs(y'*((-2*l*eye(n)-tau*exp(-tau*l)*A1)*x)));


%% Create starting vectors:
u=randn(n,1);
u=nep.Mt_solve(u);
v=randn(n,1);
v=v/(u'*nep.M_lin_comb(v,1));

m=62 % number of iterations
k1=49 % number of intermediate iterations
%% Compute reference solution (with infinite Arnoldi and Newton)
k0=100;
result_iar=iar(nep,v,k0); 
fprintf('Doing Newton corrections to have reference solutions\n');
ref_lambda=NaN*zeros(m-1,1);
ref_v_r=NaN*zeros(n,m-1);
ref_v_l=NaN*zeros(n,m-1);
[V,D,W]=eig_with_left(full(result_iar.H(1:k0,1:k0)));
vstar_arn = 1./diag(D).';

for j=1:k0
    fprintf([' ',num2str(j),':']);
    v_r=result_iar.X(:,j);
    ref_v_r(:,j)=v_r;             
    ref_lambda(j)=vstar_arn(j);                   
    ref_norms(j)=inf;
    m0=norm(nep.M(ref_lambda(j)),1);
    norm0=norm(nep.M(ref_lambda(j))*ref_v_r(:,j))/m0;
    if (norm0<1e-3)        % only do newton correction if we have a
                           % reasonable approximation
        for ii=1:5 % do some newton steps:
            fprintf('N'); % A newton step
            vv_r=ref_v_r(:,j);
            l=ref_lambda(j);
            Fp=[nep.M(l),nep.Mp(l)*vv_r; vv_r',0];
            F=[nep.M(l)*vv_r;vv_r'*vv_r-1];
            dx=-Fp\F;
            vv_r=vv_r+dx(1:n); ref_v_r(:,j)=vv_r/norm(vv_r);
            ref_lambda(j)=ref_lambda(j)+dx(end); 
            if (norm(dx)<eps*100)
                break; 
            end
        end
        norm1=norm(nep.M(ref_lambda(j))*ref_v_r(:,j))/m0;
        if (norm1>1e-8)
            % no improvement achieved with newton
            fprintf('0');
            % do not save approximation            
            ref_v_r(:,j)=NaN*v_r;             
            ref_lambda(j)=NaN*vstar_arn(j);                          
            ref_norms(j)=inf;        
        else
            ref_norms(j)=norm1;
            fprintf('+'); % Newton had positive effect
        end
    else
        % Don't do anything if initial guess is bad
        fprintf('-');
        ref_lambda(j)=NaN;        
        ref_v_r(:,j)=NaN*ones(n,1);
    end         
end
fprintf('\n');

% Remove duplicates
II=find(abs(diff(ref_lambda))<1e-8);
ref_lambda(II)=NaN*II;
ref_lambda=sortrows(ref_lambda);
ref_lambda(find(isnan(ref_lambda))) = [];

 
%% Run it
display('Running infinite bi-Lanczos');
tic
result_bilan=infbilanczos(nep,v,u,m);
% [alpha,beta,gamma,Q1,Qt1,itertime]=infbilanczos(nep,v.Y,u.Y,m);
ttoc=toc;
H_bilan=result_bilan.T; 
% H_bilan = spdiags([beta(1:m-1) alpha(1:m-1) [0;gamma(1:m-2)]], -1:1, m-1, m-1);

display('Creating convergence details eigenpairs');

%% Construct the history vector
hist_bilan=[];
for k = m-1:-1:1
    %Compute eigentriplets of matrix H.
    [V,D,W] = eig_with_left(full(H_bilan(1:k,1:k)));
    vstar = 1./diag(D).';
    
    %Sort the eigenvalues aligned to reference Ritzvalues of IAR.        
    [I,vstar]=reordervec(ref_lambda, vstar);
    hist_bilan=[hist_bilan;vstar]; %  expand hist
    V1 = V(:,find(~isnan(I))); 
    W1 = W(:,find(~isnan(I)));
    ell = size(V1,2);
    for j = 1:ell
%         v_r = Q1(:,1:k)*V1(:,j);
%         v_l = Qt1(:,1:k)*W1(:,j);
        v_r = result_bilan.Q_basis(:,1:k)*V1(:,j);
        v_l = result_bilan.Qt_basis(:,1:k)*W1(:,j);
        RES1(k,j)=nep.resnorm_lx(vstar(j),v_r);
        CON1(k,j)=nep.cond_lxy(vstar(j),v_r,v_l);
        FWE(k,j)=RES1(k,j)*CON1(k,j);
    end

    %Obtain information for the residual norms
%     for kk = 1:ell
%         VV(k,kk) = V1(k,kk);
%         WW(k,kk) = W1(k,kk);
%         Rnorm1(k,kk) = abs(V(k,kk)*beta(k));
%         Rnorm2(k,kk) = abs(W(k,kk)*gamma(k));
%     end
end
vstar_bilan = hist_bilan(1,:);
vstar_bilan_k1 = hist_bilan(m-k1,:);

errv_bilan=hist_bilan-ones(m-1,1)*ref_lambda.';
errv_bilan=flipud(errv_bilan);
% Find eigenvalues more accurate than TOL.
II1=abs(errv_bilan(end,:))<1e-3;
II_k1=abs(errv_bilan(k1,:))<1e-6;

%% Results from IAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the hist matrix: reorder eigenvalue approximations
hist_arn=[];
H_arn=result_iar.H;
for k = k0:-1:1
    [V,D,W] = eig_with_left(full(H_arn(1:k,1:k)));
    vstar = 1./diag(D).';
    % Sort the eigenvalues aligned to reference Ritzvalues of IAR.        
    [YY,II]=sortrows(abs(vstar)'); vstar=vstar(II);
    [I,vstar]=reordervec(ref_lambda, vstar);
    % expand hist
    hist_arn=[hist_arn;vstar];
end

%% compute the error 
errv_arn=hist_arn-ones(k0,1)*ref_lambda.';
errv_arn=abs(errv_arn);
errv_arn=flipud(errv_arn);

%% plot it
% vstar_k1 = 1./eig(full(H_bilan(1:k1,1:k0)));
figure(1)
plot(vstar_bilan_k1, 'b* ')
hold on;
plot(vstar_bilan_k1(II_k1), 'ro ')
title('Ritzvalues')
xlabel('Real')
ylabel('Imag')
legend('Ritz values', 'Converged Ritz values')

figure(2);
q1=semilogy(errv_arn(1:m,:),'k--');
hold on;
q2=semilogy(abs(errv_bilan(:,II1)), 'r*-');
xlabel('iteration')
ylabel('eigenvalue error')
l=legend([q2(1),q1(1)],...
         ['Infinite bi-Lanczos'],...
         ['IAR']);

figure(3);
p1=semilogy(result_iar.timing_iteration(1:m)*ones(1,size(II1,2)),errv_arn(1:m,:),'k--');
hold on;
% p2=semilogy(itertime(2:end)*ones(1,size(find(nonzeros(II1)),1)),abs(errv_bilan(:,II1)), 'r*-');
p2=semilogy(result_bilan.timing_iteration*ones(1,size(find(nonzeros(II1)),1)),abs(errv_bilan(:,II1)), 'r*-');
xlabel('wall time')
ylabel('eigenvalue error')
l=legend([p2(1), p1(1)],...
         ['Infinite bi-Lanczos'],...
         ['IAR']);


%% Print condition numbers of the approximate eigenvalues (also of
%% not converged eigenvalues).
 
[V,D,W] = eig_with_left(full(H_bilan(1:k1,1:k1)));
vstar = 1./diag(D).';
for j=1:k1
%     v_r = Q1(:,1:m-1)*V(:,j);
%     v_l = (Qt1(:,1:m-1)*W(:,j));
    v_r = result_bilan.Q_basis(:,1:k1)*V(:,j);
    v_l = (result_bilan.Qt_basis(:,1:k1)*W(:,j));
    CON(j) = nep.cond_lxy(vstar(j),v_r,v_l);
end
disp([abs(vstar)' CON']);



fprintf('Total time for IAR: %f\n',result_iar.timing_iteration(end))
fprintf('Total time for infinite bilanczos: %f\n',result_bilan.timing_iteration(end))

fprintf('Scalar product time: %f\n',sum(result_bilan.timing_scalarprod))
