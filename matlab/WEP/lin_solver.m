function [x, resvec] = lin_solver( nep, sigma, b, x0, P, m)
%LIN_SOLVER
%   Solves the system nep.M x = b
%   associated with the WEP. Uses
%   the Schur-complement approach
%   with preconditioned iterative solves


nx=nep.nx;  nz=nep.nz;
C1 = nep.C1;  C2T = nep.C2T;  Pinv = @(x) nep.Pinv(sigma,x);

b1=b(1:nx*nz);
b2=b(nx*nz+1:end);

rhs = b1-C1*Pinv(b2);
S = @(x) nep.S(sigma, x);

%tol = -Inf;
tol = 1e-12;

[x1, ~, ~, ~, resvec] = gmres(S, rhs, [], tol, m, P, [], x0); %X = GMRES(A,B,RESTART,TOL,MAXIT,M)

x2 = Pinv(b2-C2T*x1);

x=[x1; x2];


end

