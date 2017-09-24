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
                         tol=eps(real(T))*100,
                         maxit=10,
                         σ=zero(complex(T)),
                         displaylevel=0,
                         linsolvercreator::Function=backslash_linsolvercreator,
                         k=3, # Number of eigenvals to compute
                         radius=1, # integration radius
                         quad_method=:quadg_parallel, # which method to run. :quadg, :quadg_parallel, :quadgk
                         N=1000)  # Nof quadrature nodes 

    
    g=t -> radius*exp(1im*t)
    gp=t -> 1im*radius*exp(1im*t)

    n=size(nep,1);
    srand(10); # Reproducability
    Vh=randn(n,k)
    if (k>n)
        println("k=",k," n=",n);
        error("Cannot compute more eigenvalues than the size of the NEP with contour_beyn()");

    end

    function local_linsolve(λ,V)
        @ifd(print("."))
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        return lin_solve(M0inv,V);
    end

    # Constructing integrands
    Tv0= λ ->  local_linsolve(λ,Vh)
    Tv1= λ -> λ*Tv0(λ)
    f1= t-> Tv0(g(t))*gp(t) 
    f2= t -> Tv1(g(t))*gp(t)
    @ifd(print("Computing integrals"))

    # Naive version, where we compute two separate integrals

    local A0,A1
    if (quad_method == :quadg_parallel)
        @ifd(print(" using quadg_parallel"))
        A0=quadg_parallel(f1,0,2*pi,N);
        A1=quadg_parallel(f2,0,2*pi,N);        
    elseif (quad_method == :quadg)
        @ifd(print(" using quadg"))
        A0=quadg(f1,0,2*pi,N);
        A1=quadg(f2,0,2*pi,N);
    elseif (quad_method == :quadgk)
        @ifd(print(" using quadgk"))
        A0,tmp=quadgk(f1,0,2*pi,reltol=tol);
        A1,tmp=quadgk(f2,0,2*pi,reltol=tol);
    else
        error("Unknown quadrature method:"*String(quad_method));
    end
    @ifd(println("."));
    # Don't forget scaling
    A0=A0/(2im*pi);
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

#  Carries out Gauss quadrature (with N) discretization points
#  by call to @parallel
function quadg_parallel(f,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    S=@parallel (+) for i = 1:N
        f(t[i])*w[i]
    end
    return S;
end


function quadg(f,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    S=zeros(size(f(t[1])))
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    for i = 1:N
        S+= f(t[i])*w[i]
    end
    return S;
end




