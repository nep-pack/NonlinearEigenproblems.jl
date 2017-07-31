function [x,rho,res] = lopcg(A,x,max_it)
%RQMIN1 [x,rho] = lopcg(A,M,x0,tol,C)
% Locally Optimal Proconditioned CG algorithm for
% computing the smallest eigenvalue of A*x = lambda*M*x,f
% where A and M are symmetrisch, M spd.
% x0 initial vektor
% C’*C preconditioner
% tol: stopping criterion:
% (C’*C)\(A*x - lam*M*x) < tol
% PA 2002-07-3
n =size(A,1);
AA=A'*A;
u =x;
q =sqrt(x'*u);
x =x/q; u = u/q;
v =AA*x;
rho= x'*v;

p = zeros(n,0); 



k = 1;
while k < max_it
    g = v - rho*u;
    % gradient
    res(k,1)=norm(A*x);
    


    aa = [x -g p]'*[v -AA*g AA*p];

    
    aa = (aa+aa')/2;
    mm = [x -g p]'*[u -g p]; 
    mm = (mm+mm')/2;
    [qq,ll] = eig(aa,mm);
    [rho,ii] = min(diag(ll));
    delta = qq(:,ii);
    p = [-g p]*delta(2:end);
    x = delta(1)*x + p;
    u = x;
    q = sqrt(x'*u);
    x = x/q; u = u/q;
    v = AA*x;
    k = k + 1;
    
end

end
