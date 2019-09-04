using NonlinearEigenproblems, Random, LinearAlgebra


nep=nep_gallery("dep1")
# iar linearization
d=5
M=diagm( 0 =>  ones(d) )[2:end,:]
N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]

Av=Array{AbstractMatrix,1}(undef, d)
Bv=Array{AbstractMatrix,1}(undef, d)
for j=1:d
    Av[j]=compute_Mder(nep,0,j)/j
end
Bv[1]=-compute_Mder(nep,0,0)
for j=2:d
    Bv[j]=zero(Av[1])
end

cp=CORK_pencil(M,N,Av,Bv)
AA,BB=build_CORK_pencil(cp)
