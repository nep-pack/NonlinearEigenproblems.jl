using NonlinearEigenproblems, Random, LinearAlgebra, PyPlot, Revise
#include("extract_nleigs_coefficients.jl")
#include("nleigs_coefficients.jl")

nep=nep_gallery("dep0");
unit_square = float([1+1im, 1-1im, -1-1im,-1+1im])

r=5; θ=range(0,stop=2π,length=1000); Σ=r*cos.(θ) + 1im*r*sin.(θ); Ξ=[-10.0]
# nleigs linearization
cp=compute_CORK_pencil(nep,NleigsCorkLinearization(Σ=Σ,Ξ=[-10.0]))

AA,BB=build_CORK_pencil(cp)


λ2=eigvals(Matrix(AA),Matrix(BB))

λ,_=iar(nep,maxit=200,tol=1e-8,neigs=Inf)
pygui(true)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:red,linestyle=:none)
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
