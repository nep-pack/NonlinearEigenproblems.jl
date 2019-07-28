# Evaluate a ChebPEP in a point: <=> Interpolate
# Compute coefficients:
# http://inis.jinr.ru/sl/M_Mathematics/MRef_References/Mason,%20Hanscomb.%20Chebyshev%20polynomials%20(2003)/C0355-Ch06.pdf


f=s-> ones(3)-(1:3)*sin(s);
#f=s-> 1-1/(s+6);
#f=s-> 1+s;
k=3;
a=-1;b=4;

(Fk,xk)=chebyshev_eval(a,b,k,f);

Ck0=chebyshev_compute_coefficients(a,b,Fk,xk);
#Ck1=chebyshev_compute_coefficients_naive(a,b,Fk,k);

#Ck=Ck1;
xk=get_chebyshev_nodes(Float64,a,b,k)
xx=xk[1];
ff0=mapreduce(i->Ck0[i]*cheb_f(a,b,xx,i-1), +, 1:k)
#ff1=mapreduce(i->Ck1[i]*cheb_f(a,b,xx,i-1), +, 1:k)

@show norm(ff0-f(xx))
#@show ff1-f(xx)




nep=nep_gallery("dep0",5);
nep2=ChebPEP(nep,11,-1,1);
nep3=nep2.spmf;

s=0;

compute_Mder(nep,s)-compute_Mder(nep3,s)

λ1=iar(nep,neigs=2,logger=1)[1]
λ2=iar(nep2,neigs=2,logger=1,errmeasure=ResidualErrmeasure)[1]
λ3=iar(nep3,neigs=2,logger=1,errmeasure=ResidualErrmeasure)[1]


@show norm(λ1-λ3)





# Implement linearization technique (e.g. Effenberger & Kressner)
#function polyeig(pep::ChebPEP)
#
#end


(λ,V)=polyeig(nep2)

norm(compute_Mder(nep2,λ[1])*V[:,1])
