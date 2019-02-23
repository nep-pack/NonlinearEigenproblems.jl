# Beyn contour integral approach

#using QuadGK # Disabled to lessen requirements
using Distributed
using LinearAlgebra
using Random

export contour_beyn

"""
    λv,V=contour_beyn([eltype],nep;[tol,][displaylevel,][σ,],[linsolvercreator,][k,][radius,][quad_method,][N,][neigs,])

The function computes eigenvalues using Beyn's contour integral approach,
using a circle centered at `σ` with radius `radius`. The quadrature method
is specified in `quad_method` (`:ptrapz`, `:quadg`,`:quadg_parallel`,`:quadgk`). `k`
specifies the number of computed eigenvalues. `N` corresponds to the
number of quadrature points. Circles are the only supported contours. The
`linsolvercreator` must create a linsolver that can handle (rectangular) matrices
as right-hand sides, not only vectors.

The kwargs `neigs` specifies the number of wanted eigvals, and `k` is the number
of columns in the subspace (default `k=neigs+1`). If you give the `k` parameter and set `neigs=typemax(Int)` all found eigenvalues will be returned. The kwarg `sanity_check`
decides if checking of eigpars should be done. If disabled, the method
returns `k` (potentially inaccurate) eigpairs. The parameters `errmeasure` and
`tol` are used for the sanity check.

# Example
```julia-repl
julia> using LinearAlgebra
julia> nep=nep_gallery("dep0");
julia> λv,V=contour_beyn(nep,radius=1,neigs=2,quad_method=:ptrapz);
julia> norm(compute_Mlincomb(nep,λv[1],V[:,1])) # Eigenpair 1
5.778617503485546e-15
julia> norm(compute_Mlincomb(nep,λv[2],V[:,2])) # Eigenpair 2
3.095638020248726e-14
```
# References
* Wolf-Jürgen Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications 436 (2012) 3839–3863
"""
contour_beyn(nep::NEP;params...)=contour_beyn(ComplexF64,nep;params...)
function contour_beyn(::Type{T},
                         nep::NEP;
                         tol::Real=sqrt(eps(real(T))), # Note tol is quite high for this method
                         σ::Number=zero(complex(T)),
                         displaylevel::Integer=0,
                         linsolvercreator::Function=backslash_linsolvercreator,
                         neigs::Integer=2, # Number of wanted eigvals
                         k::Integer=neigs+1, # Columns in matrix to integrate
                         radius::Real=1, # integration radius
                         quad_method::Symbol=:ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz
                         N::Integer=1000,  # Nof quadrature nodes
                         errmeasure::Function =
                           default_errmeasure(nep::NEP),
                         sanity_check=true
                        )where{T<:Number}

    # Geometry
    g=t -> radius*exp(1im*t)
    gp=t -> 1im*radius*exp(1im*t) # Derivative

    n=size(nep,1);

    if (k>n)
        error("Cannot compute more eigenvalues than the size of the NEP with contour_beyn() k=",k," n=",n);

    end
    if (k<=0)
        error("k must be positive, k=",k,
              neigs==typemax(Int) ? ". The kwarg k must be set if you use neigs=typemax" : ".")
    end


    Random.seed!(10); # Reproducability
    Vh=Array{T,2}(randn(real(T),n,k)) # randn only works for real




    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
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
        @ifd(print(" using quadg_parallel"))
        A0=quadg_parallel(f1,0,2*pi,N);
        A1=quadg_parallel(f2,0,2*pi,N);
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
    A0[:,:] = A0 ./(2im*pi);
    A1[:,:] = A1 ./(2im*pi);

    @ifd(print("Computing SVD prepare for eigenvalue extraction "))
    V,S,W = svd(A0)
    V0 = V[:,1:k]
    W0 = W[:,1:k]
    B = (copy(V0')*A1*W0) * Diagonal(1 ./ S[1:k])

    rank_drop_tol=tol;
    p = count( S/S[1] .> rank_drop_tol);

    @ifd(println(" p=",p));

    # Extract eigenval and eigvec approximations according to
    # step 6 on page 3849 in the reference
    @ifd(println("Computing eigenvalues "))
    λ,VB=eigen(B)
    λ[:] = λ .+ σ

    @ifd(println("Computing eigenvectors "))
    V = V0 * VB;
    for i = 1:k
        normalize!(V[:,i]);
    end

    if (!sanity_check)
        return (λ,V)
    end

    # Compute all the errors
    errmeasures=zeros(real(T),k);
    for i = 1:k
        errmeasures[i]=errmeasure(λ[i],V[:,i]);
    end

    good_index=findall(errmeasures .< tol);

    # Index vector for sorted to distance from σ
    sorted_good_index=
       good_index[sortperm(map(x->abs(σ-x), λ[good_index]))];

    # Remove all eigpairs not sufficiently accurate
    # and potentially eigenvalues we do not want.
    local Vgood,λgood
    if( size(sorted_good_index,1) > neigs)
        @ifd(println("Removing unwanted eigvals: neigs=",neigs,"<",size(sorted_good_index,1),"=found_eigvals"))
        Vgood=V[:,sorted_good_index[1:neigs]];
        λgood=λ[sorted_good_index[1:neigs]];
    else
        Vgood=V[:,sorted_good_index];
        λgood=λ[sorted_good_index];
    end


    if (p==k)
       @warn "Rank-drop not detected, your eigvals may be correct, but the algorithm cannot verify. Try to increase k." S
    end

    if (size(λgood,1)<neigs  && neigs < typemax(Int))
       @warn "We found less eigvals than requested. Try increasing domain, or decreasing `tol`." S
    end

    return (λgood,Vgood)
end

#  Carries out Gauss quadrature (with N) discretization points
#  by call to @parallel
function quadg_parallel(f,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    S=@distributed (+) for i = 1:N
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
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S = zero(f(t[1]))
    for i=1:N
        S+=f(t[i])
    end
    return h*S;
end
