function [eigval, eigvec] = main_func_WEP(nz, N, lambda, delta, waveguide_str)
% RESIDUAL INVERSE ITERATION FOR WEP
% IMPLEMENTING INVERSE WITH SHCUR COMPLEMENT
% SCHUR-COMPLEMENT IS SOLVED WITH GMRES
%
% Sylvester-based preconditioning for the waveguide eigenvalue problem;
% Emil Ringh, Giampaolo Mele, Johan Karlsson, and Elias Jarlebring
%
% NOTATION:
%      lambda = initial guess of the eigenvalue
%      v = initial guess of the eigenvector
%      sigma = fixed shift used in Resinv
%      mu = the new iterate of the eigenvalue
%      w = new iterate of the eigenvector


% number of discretization points
%n = 105;
n = nz;

% Number of "interior" blocks   !!!OBS: Must have n/N = integer
%N = 35;

% number of max iterations in GMRES (in linear solve)
m = 100;

% number of steps in residual inverse iteration
max_nrof_resinv_steps = 30;

% backward error stopping criterion in residual inverse iteration
nep_stop_crit = 1e-14;

% INITIAL GUESSES AND SHIFT
%lambda = -0.51-0.38i;



% generate the problem
%nz = n;
nx = n + 4; %OBS: Should be nz + 4.

nep_options.delta = delta;
nep_options.wg = waveguide_str;    %  CHALLENGE  TAUSCH
nep = nep_wg_generator(nx, nz, nep_options);


sigma = lambda;
mu = NaN(max_nrof_resinv_steps+1,1);
mu(1) = lambda;
seed = 1052422089; %Random seed gotten from hashing the clock at an instant
rng(seed, 'twister');
v = ones(nx*nz+2*nz,1);%rand(nx*nz+2*nz,1);
w=v/norm(v);


tic
% CONSTRUCT THE PRECONDITIONER
kk = mean(nep.K(:));
K = nep.K - kk;

Linv=@(X) fft_wg( X, sigma, kk, nep.hx, nep.hz );
Pm_inv=@(x) -nep.Pm_inv(sigma, x);              %OBS! The minus sign!
Pp_inv=@(x) -nep.Pp_inv(sigma, x);              %OBS! The minus sign!
dd1 = nep.d1/nep.hx^2;
dd2 = nep.d2/nep.hx^2;
[ M ] = generate_smw_matrix( n, N, Linv, dd1, dd2, Pm_inv, Pp_inv, K, true );
precond =@(x) reshape(solve_smw( M, reshape(x,nz,nx), Linv, dd1, dd2, Pm_inv, Pp_inv, K ),nz*nx,1);
t_precalc = toc;


% RESIDUAL INVERSE ITERATION
nep_backward_error = NaN(max_nrof_resinv_steps+1, 1);
nep_resnorm  = NaN(max_nrof_resinv_steps+1, 1);
linsys_gmres_res = cell(1,max_nrof_resinv_steps);
linsys_resnorm  = NaN(max_nrof_resinv_steps+1, 1);

nep_resnorm(1) = norm(nep.M(mu(1),w));
fprintf('NEP Residual norm at start = %d\n', nep_resnorm(1));
nep_backward_error(1) = backward_error(nep, mu(1), w);
fprintf('NEP backward error at start = %d\n', nep_backward_error(1));
dw = zeros(nx*nz+2*nz, 1);
idx = 0;
t_resinv = NaN(max_nrof_resinv_steps, 1);


while((idx<max_nrof_resinv_steps) && nep_backward_error(idx+1)>nep_stop_crit)
tic
    idx = idx +1;
    fprintf('%i out of maximum %i\n', idx, max_nrof_resinv_steps);

    % 1) New iterate of the eigenvalue solving Rayleigh quotent
    f  = @(gamma) w'*nep.M(gamma,w);
    fp = @(gamma) w'*nep.p_M(gamma,w);
    mu(idx+1) = newton(f, fp, mu(idx));
    
    
    % 2) Compute the residual 
    r = nep.M(mu(idx+1),w);
    
    
    % 3) Compute correction to the eigenvector
    x0 = zeros(nx*nz,1);
    [dw, resvec] = lin_solver(nep, sigma, r, x0, precond, m);
    fprintf('    GMRES preconditioned final residual norm = %d\n', resvec(end));
    linsys_resnorm(idx +1) = norm(nep.M(sigma,dw) - r);
    fprintf('    Lin-System Actual residual norm for M = %d\n', linsys_resnorm(idx +1));
    linsys_gmres_res{idx} = resvec;
    
    
    % 4) Update eigenvector approximation
    w = w-dw;
    w = w/norm(w);
    
    
    % 5) Compute the backward error
    nep_resnorm(idx+1) = norm(nep.M(mu(idx+1),w));
    fprintf('    NEP Residual norm = %d\n', nep_resnorm(idx+1))
    nep_backward_error(idx+1) = backward_error(nep, mu(idx+1), w);
    fprintf('    NEP backward error = %d\n', nep_backward_error(idx+1));
    
t_resinv(idx) = toc;
end

eigval = mu(idx+1);
eigvec = w;
end
%fprintf('Computation time in seconds: %s\n', num2str(sum(t_resinv(~isnan(t_resinv)))+t_precalc))
%
%
%figure
%for idx = 1:max_nrof_resinv_steps
%    semilogy(linsys_gmres_res{idx},'--*')
%    hold on
%end
%xlabel('GMRES iterations')
%ylabel('Preconditioned relative residual norm')
%
%figure
%semilogy(0:max_nrof_resinv_steps, nep_backward_error, '-*')
%xlabel('Number of steps in Residual inverse iteration')
%ylabel('Backwards error estimate')
%hold on
%dd = 0:(find(~isnan(nep_backward_error),1,'last')-1);
%dist = abs(mu(find(~isnan(mu),1,'last'))-sigma);
%lin_dec = 1e-3*dist.^dd;
%semilogy(dd,lin_dec)
%legend('Experiment', 'Prediction')

%W = [w(nx*nz+1:nx*nz+nz), reshape(w(1:nx*nz),nz,nx), w(nx*nz+nz+1:nx*nz+2*nz)]/max(max(abs(w)));
%x = linspace(nep.xm, nep.xp, nep.nx);
%z = linspace(0, 1 , nep.nz);
%figure
%imagesc(x, z, abs(W)); colorbar
%set(gca,'YDir','normal')
%daspect([1 1 1])
