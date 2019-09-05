using NonlinearEigenproblems, Random, LinearAlgebra, PyPlot, Revise
include("extract_nleigs_structure.jl")


nep=nep_gallery("dep0");
unit_square = float([1+1im, 1-1im, -1-1im,-1+1im])

r=5
θ=range(0,stop=2π,length=1000)
Σ=r*cos.(θ) + 1im*r*sin.(θ)
D, β, ξ, σ=nleigs_structure(nep,unit_square);

d=19
σ=σ[1:d+1]; β=β[1:d+1]; ξ=ξ[1:d+1]

M=diagm( -1 => σ[1:d], 0 =>  β[1:d] )[2:end-1,1:end-1]
Av=Array{AbstractMatrix,1}(undef, d)
Av[1:d-1]=D[1:d-1]; Av[d]=D[d]-σ[d]/β[d+1]*D[d+1]


N=diagm( -1 => ones(d), 0 =>  β[1:d]./ξ[1:d] )[2:end-1,1:end-1]
Bv=Array{AbstractMatrix,1}(undef, d)
Bv[1:d-1]=D[1:d-1]/ξ[d+1]; Bv[d]=D[d]/ξ[d+1]-D[d+1]/β[d+1]
ck=CORK_pencil(M,N,Av,Bv)

AA,BB=build_CORK_pencil(ck)

λ2=eigvals(Matrix(AA),Matrix(BB))

#λ,_=iar(nep,maxit=200,tol=1e-8,neigs=Inf)
pygui(true)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:red,linestyle=:none)
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
