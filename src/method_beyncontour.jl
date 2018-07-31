# Beyn contour integral approach

#using QuadGK # Disabled to lessen requirements 

export contour_beyn

"""
    λv,V=contour_beyn([eltype],nep;[errmeasure,][tol,][maxit,][displaylevel,][σ,],[linsolvercreator,][k,][radius,][quad_method,][N])

The function computes eigenvalues using Beyn's contour integral approach,
using a circle centered at `σ` with radius `radius`. The quadrature method
is specified in `quad_method` (`:ptrapz`, `:quadg`,`:quadg_parallel`,`:quadgk`). `k`
specifies the number of computed eigenvalues. `N` corresponds to the
number of quadrature points. 

# Example
```julia-repl 
julia> nep=nep_gallery("dep0");
julia> λv,V=contour_beyn(nep,radius=1,k=2,quad_method=:ptrapz);
julia> minimum(svdvals(compute_Mder(nep,λv[1])))
1.6106898471314257e-16
```
# References
* Wolf-Jürgen Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications 436 (2012) 3839–3863
"""


contour_beyn(nep::NEP;params...)=contour_beyn(Complex128,nep;params...)
function contour_beyn(::Type{T},
                         nep::NEP;
                         errmeasure::Function =
                         default_errmeasure(nep::NEP),
                         tol::Real=eps(real(T))*100,
                         maxit::Integer=10,
                         σ::Number=zero(complex(T)),
                         displaylevel::Integer=0,
                         linsolvercreator::Function=backslash_linsolvercreator,
                         k::Integer=3, # Number of eigenvals to compute
                         radius::Real=1, # integration radius
                         quad_method::Symbol=:ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz
                         N::Integer=1000  # Nof quadrature nodes 
                         )where{T<:Number}
    
    g=t -> radius*exp(1im*t)
    gp=t -> 1im*radius*exp(1im*t)

    n=size(nep,1);
    srand(10); # Reproducability
    Vh=Array{T,2}(randn(real(T),n,k)) # randn only works for real
    
    if (k>n)
        println("k=",k," n=",n);
        error("Cannot compute more eigenvalues than the size of the NEP with contour_beyn()");

    end

    function local_linsolve(λ::TT,V::Array{TT,2}) where {TT<:Number}
        @ifd(print("."))
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end

    # Constructing integrands
    Tv0= λ ->  local_linsolve(T(λ),Vh)
    Tv1= λ -> λ*Tv0(λ)
    f1= t-> Tv0(g(t))*gp(t) 
    f2= t -> Tv1(g(t))*gp(t)
    @ifd(print("Computing integrals"))

    # Naive version, where we compute two separate integrals

    local A0,A1
    if (quad_method == :quadg_parallel)
        #@ifd(print(" using quadg_parallel"))
        #A0=quadg_parallel(f1,0,2*pi,N);
        #A1=quadg_parallel(f2,0,2*pi,N);        
        error("disabled");
    elseif (quad_method == :quadg)
        #@ifd(print(" using quadg"))
        #A0=quadg(f1,0,2*pi,N);
        #A1=quadg(f2,0,2*pi,N);
        error("disabled");        
    elseif (quad_method == :ptrapz)
        @ifd(print(" using ptrapz"))
        A0=ptrapz(f1,0,2*pi,N);
        A1=ptrapz(f2,0,2*pi,N);
    elseif (quad_method == :quadgk)
        #@ifd(print(" using quadgk"))
        #A0,tmp=quadgk(f1,0,2*pi,reltol=tol);
        #A1,tmp=quadgk(f2,0,2*pi,reltol=tol);
        error("disabled");
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


# Trapezoidal rule for a periodic function f
function ptrapz(f,a,b,N)
    h=(b-a)/N
    t=linspace(a,b-h,N);
    S=zeros(f(t[1]))
    for i=1:N
        S+=f(t[i])
    end
    return h*S;
end




