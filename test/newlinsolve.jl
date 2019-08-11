using SparseArrays;



n=5;
α=0.1;
A=spdiagm(0=>ones(n),1=>α*ones(n-1),-1=>α*ones(n-1));
B=spdiagm(0=>ones(n));
C=spdiagm(0=>(1:n)/n);

nep= SPMF_NEP([A,B,C],[s->one(s),s->s,s->exp(s)])

quasinewton(nep,λ=-1.0,v=ones(n),logger=1);



quasinewton(nep,λ=-1,v=ones(n),logger=1,linsolvercreator=BackslashLinSolverCreator());


creator=GMRESLinSolverCreator(log = true,verbose=true)
quasinewton(nep,λ=-1,v=ones(n),logger=1,
            linsolvercreator=creator);
