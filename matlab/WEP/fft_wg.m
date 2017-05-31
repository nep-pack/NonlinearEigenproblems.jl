function [ X ] = fft_wg( C, gamma, kk, hx, hz )
%FFT_WG SOLVE THE SYLVESTER EQUATION ASSOCIATED TO THE WEP
%
%   SOLVE THE SYLVESTER MATRIX EQUATION
%   (Dzz+2*gamma*Dz+(gamma^2+kk)*I)*U + U*Dxx = C
%   we define
%   A=Dzz+2*gamma*Dz+(gamma^2+kk)*I
%   B=Dxx
%   X=U
%   and we solve the Sylvester equation
%   A*X+X*B=C
%
%   EXAMPLE OF EXECUTION
% % clear; close all; clc;
% % nz = 100;  nx = nz+4;  hx=5/nx;  hz=2/nz;
% % % BUILD THE DISCRETIZATION MATRICES
% % ex = ones(nx,1);  Dxx = spdiags([ex -2*ex ex], -1:1, nx, nx);
% % ez = ones(nz,1);  Dzz = spdiags([ez -2*ez ez], -1:1, nz, nz);
% % Dz  = spdiags([-ez ez], [-1 1], nz, nz);
% % 
% % % IMPOSE PERIODICITY
% % Dz(1,end) = -1;     Dz(end,1) = 1;  Dzz(1, end) = 1;    Dzz(end, 1) = 1;
% % 
% % % SCALE THE MATRICES
% % Dxx = Dxx/(hx^2);   Dzz = Dzz/(hz^2);   Dz  = Dz/(2*hz);
% % 
% % % BUILD THE SYLVESTER EQUATION
% % C=rand(nz,nx);  kk=rand;    gamma=-rand-rand*1i;    A=Dzz+(2*gamma)*Dz+(gamma^2+kk)*speye(nz,nz);  B=Dxx; 
% % 
% % % SOLVE THE SYLVESTER EQUATION
% % tt=cputime;     X = fft_wg( C, gamma, kk, hx, hz ); tt=cputime-tt;
% % fprintf('time fft_wg %f',tt);
% % 
% % % SOLVE THE SYLVESTER EQUATION WITH LYAP
% % tt=cputime;    XX=lyap(A,B,-C); tt=cputime-tt;
% % fprintf('\ntime lyap %f ',tt);
% % 
% % fprintf('\nresidual fft_wg %d',norm(A*X + X*B - C));
% % fprintf('\nresidual lyap %d',norm(A*XX + XX*B - C));
% % fprintf('\ndifference fft_wg and lyap %d\n',norm(X-XX));


nz = size(C,1);
nx = size(C,2);

alpha = gamma^2+kk;

% eigenvalues of A
v=zeros(nz,1);   v(1)=-2;    v(2)=1;     v(end)=1;   v=v/(hz^2);
w=zeros(nz,1);   w(2)=1;     w(end)=-1;   w=w*(gamma/hz);
D=fft(v+w)+alpha;

V=@(X)  fft(X)/sqrt(nx);
Vh=@(X) ifft(X)*sqrt(nx);

% eigenvalues of B
S=-(4/hx^2)*sin(pi*(1:nx)/(2*(nx+1))).^2;  S=S.';

% solve the diagonal matrix equation
Z=zeros(nz,nx);

CC=Vh( Wh(C')' );

for k=1:nx
    Z(:,k)=CC(:,k)./(D+S(k));
end


% change variables
X=V((W(Z'))');

end

function [ WX ] = W( X )
%W Compute the action of the matrix W
%   W is the matrix of the eigenvectors of the second derivative Dxx
%   W*X can be computed with FFTs

WX=(1i/2)*(F(X)-Fh(X));

nz = size(X,1);    WX=WX/sqrt((nz+1)/2);


end

function [ WX ] = Wh( X )
%Wh Compute the action of the matrix Wh
%   Wh is the transpose of the matrix of the eigenvectors of the second derivative Dxx
%   Wh*X can be computed with FFTs


WX=(Fh(X)-F(X))/(2i);

nz = size(X,1);    WX=WX/sqrt((nz+1)/2);


end

function [ v ] = F( v )
%F is an auxiliary function
%   Compute the FFT of a matrix

m=size(v,2);

v=[zeros(1,m); v];
n=size(v,1);

v=fft(v, 2*n);
v=v(2:n,:);

end

function [ v ] = Fh( v )
%Fh is an auxiliary function
%   Compute the IFFT of a matrix

m=size(v,2);

v=[zeros(1,m); v];
n=length(v);

v=ifft(v, 2*n);
v=v(2:n,:);
v=v*2*n;

end
