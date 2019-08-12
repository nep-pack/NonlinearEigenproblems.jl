using SparseArrays,LinearAlgebra;



n=100;
α=0.01;
A=spdiagm(0=>ones(n),1=>α*ones(n-1),-1=>α*ones(n-1));
B=spdiagm(0=>ones(n));
C=spdiagm(0=>(1:n)/n);

nep= SPMF_NEP([A,B,C],[s->one(s),s->s,s->exp(s)])

λ0=-1.002
quasinewton(nep,λ=λ0,v=ones(n),logger=1);


#quasinewton(nep,λ=λ0,v=ones(n),logger=1,linsolvercreator=BackslashLinSolverCreator());


D0=(Diagonal(compute_Mder(nep,λ0)));
creator=GMRESLinSolverCreator(Pl=D0, tol=1e-1, log=true)
(λ,x)=quasinewton(nep,λ=λ0,v=ones(n),logger=1,
                  linsolvercreator=creator,tol=1e-16,maxit=100);
normalize!(x)
norm(compute_Mlincomb(nep,λ,x))
