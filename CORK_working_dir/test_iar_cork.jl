using NonlinearEigenproblems, Random, LinearAlgebra, PyPlot, SparseArrays

pygui(true)

#nep=nep_gallery("dep0_sparse",100,0.1)
n=10; A0=rand(n,n); A0=-A0'*A0; A1=rand(n,n)
nep=DEP([A0,A1],[0,1.0])

# iar linearization
cp=compute_CORK_pencil(nep,IarCorkLinearization(d=40))

AA,BB=build_CORK_pencil(cp)

λ2=eigvals(Matrix(AA),Matrix(BB))

λ,_=iar(nep,maxit=200,tol=1e-8,neigs=Inf)

plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:red,linestyle=:none)
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
