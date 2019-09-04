using NonlinearEigenproblems, Random, LinearAlgebra, PyPlot


nep=nep_gallery("dep0")
# iar linearization
d=5
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

λ=eigvals(AA,BB)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)
