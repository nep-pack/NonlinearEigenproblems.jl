    clear all
close all
clc


% TEST MATRICES 1
n=200; p=5/n; 
nep.A1=sprand(n,n,p); nep.A2=sprand(n,n,p);
%e = ones(n,1);  nep.A1 = spdiags([e -2*e e], -1:1, n, n);   nep.A2=spdiags([e 2*e e], -1:1, n, n);

% precomputation (lu-factorization of M0)
[nep.L, nep.U, nep.P, nep.Q, nep.R] = lu(nep.A1+nep.A2);    nep.R=diag(nep.R);




kmax=20;   p=10;


S=rand;
Y=rand(n,1);
N_restart=5;


[ ERR_restarting ] = IAR_restarting( nep, S, Y, kmax, p, N_restart );


% PLOT THE ERROR HISTORY


close all
ERR=NaN(kmax,kmax-p,N_restart-1);
for N=1:N_restart-1    
    ERR(:,:,N)=ERR_restarting(:,p+1:kmax,N+1);
end

clear E EE

for k=1:kmax      
    E(:,k)=ERR_restarting(k,1:kmax,1);
end

for k=1:kmax    
    EE(:,k)=ERR(k,:).';
end

E(end,:)=E(end-1,:);
T=[E;EE];

for k=1:kmax
    semilogy(T(:,k),'-k');
    hold on
end

line([kmax kmax], get(gca, 'ylim'));
for N=1:N_restart-1    
     line([kmax+N*(kmax-p) kmax+N*(kmax-p)], get(gca, 'ylim'));    
end

    