using NonlinearEigenproblems, Random, LinearAlgebra, PyPlot

pygui(true)

#nep=nep_gallery("dep0_sparse",100,0.1)
n=10; A0=rand(n,n); A0=-A0'*A0; A1=rand(n,n)
nep=DEP([A0,A1],[0,1.0])


λ,_=iar(nep,maxit=200,tol=1e-8,neigs=Inf)

# iar linearization
d=10
M=diagm( 0 =>  ones(d) )[2:end,:]
N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]

Av=Array{AbstractMatrix,1}(undef, d)
Bv=Array{AbstractMatrix,1}(undef, d)

Av[1]=-compute_Mder(nep,0,0)
for j=2:d
    Av[j]=zero(Av[1])
end

for j=1:d
    Bv[j]=compute_Mder(nep,0,j)/j
end


cp=CORK_pencil(M,N,Av,Bv)
AA,BB=build_CORK_pencil(cp)

λ2=eigvals(Matrix(AA),Matrix(BB))

plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:red,linestyle=:none)
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
