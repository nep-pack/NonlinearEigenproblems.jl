using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot

#nep=nep_gallery("dep_symm_double",10)
n=20
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = pi, length = n)
h=x[2]-x[1];
#h=pi
LL=LL/(h^2)
LL=-kron(LL,LL)
A=LL

b=broadcast((x,y)->sin(x+y),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1.0])

# COMPUTE REFERENCE EIGENVALUES WITH IAR
λ,v=iar(nep;maxit=100,tol=1e-12,neigs=Inf,logger=1)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)

#Σ = float([0-3im, -3-3im, -3+3im,0+3im])
θ=range(0,stop=2π,length=1000); r=4;
Σ=r*cos.(θ) + 1im*r*sin.(θ)


λ2,v2=nleigs(nep,Σ,logger=1,tol=1e-6,minit=200,nodes=[0.0+0*im])
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
#plot(real(Σ),imag(Σ),marker="s",markerfacecolor=:none,c=:black,linestyle=:none)




#λ3,v=contour_beyn(nep,tol=1e-4,neigs=20,logger=1,N=10000,radius=3,sanity_check=false);
#plot(real(λ3),imag(λ3),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
