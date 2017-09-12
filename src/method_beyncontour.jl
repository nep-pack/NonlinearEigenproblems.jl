# Beyn contour integral approach
using QuadGK

export contour_beyn

"""
    contour_beyn computes eigenvalues using Beyn's contour integral approach
"""


contour_beyn(nep::NEP;params...)=contour_beyn(Complex128,nep;params...)
function contour_beyn{T}(::Type{T},
                         nep::NEP;
                         errmeasure::Function =
                         default_errmeasure(nep::NEP),
                         tolerance=eps(real(T))*100,
                         maxit=10,
                         σ=zero(complex(T)),
                         displaylevel=0,
                         linsolvercreator::Function=default_linsolvercreator,
                         k=3, # Number of eigenvals to compute
                         radius=1) # integration radius

    
    g=t -> radius*exp(1im*t)
    gp=t -> 1im*radius*exp(1im*t)

    n=size(nep,1);
    Vh=randn(n,k)
    if (k>n)
        println("k=",k," n=",n);
        error("Cannot compute more eigenvalues than the size of the NEP with contour_beyn()");

    end
    
    function local_linsolve(λ,V)
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        return lin_solve(M0inv,V);
    end

    # Constructing integrands
    Tv0= λ ->  local_linsolve(λ,Vh)
    Tv1= λ -> λ*Tv0(λ)
    f1= t-> Tv0(g(t))*gp(t) 
    f2= t -> Tv1(g(t))*gp(t)
    @ifd(println("Computing integrals"))

    # Naive version, where we compute two separate integrals
    A0,tmp=quadgk(f1,0,2*pi,reltol=tolerance);
    A0=A0/(2im*pi);
    A1,tmp=quadgk(f2,0,2*pi,reltol=tolerance);
    A1=A1/(2im*pi);    
    
    @ifd(println("Computing SVD prepare for eigenvalue extraction "))
    V,S,W = svd(A0);    
    V0=V[:,1:k];
    W0=W[:,1:k];
    B=(V0'*A1*W0)/diagm(S[1:k]);
    if ((maximum(S)/minimum(S))>1/sqrt(eps()))
        warn("Rank drop detected in A0. The disc probably has fewer eigenvalues than those in the disc. Try decreasing k in contour integral solver")
        println(S)
    end
    

    @ifd(println("Computing eigenvalues "))
    λ,v=eig(B);  # Eigenvector extraction not implemented yet

    return (λ+σ,NaN*v)

end





