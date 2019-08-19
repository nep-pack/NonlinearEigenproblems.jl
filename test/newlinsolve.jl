using SparseArrays,LinearAlgebra;



n=100;
α=0.01;
A=spdiagm(0=>ones(n),1=>α*ones(n-1),-1=>α*ones(n-1));
B=spdiagm(0=>ones(n));
C=spdiagm(0=>(1:n)/n);

nep= SPMF_NEP([A,B,C],[s->one(s),s->s,s->exp(s)]);

λ0=-1.02
quasinewton(nep,λ=λ0,v=ones(n),logger=1);


#quasinewton(nep,λ=λ0,v=ones(n),logger=1,linsolvercreator=BackslashLinSolverCreator());


D0=(Diagonal(compute_Mder(nep,λ0)));
creator=GMRESLinSolverCreator(Pl=D0, tol=1e-6, log=true)
(λ,x)=augnewton(nep,λ=λ0,v=ones(n),logger=1,
                  linsolvercreator=creator,tol=1e-13,maxit=1000);
@show λ
normalize!(x)
@show norm(compute_Mlincomb(nep,λ,x))

(λ,x)=augnewton(nep,λ=λ0,v=ones(n),logger=1,
                tol=1e-16,maxit=100);
@show λ
normalize!(x)
@show norm(compute_Mlincomb(nep,λ,x))
