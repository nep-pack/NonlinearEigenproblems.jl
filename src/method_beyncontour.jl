# Beyn contour integral approach

#using QuadGK # Disabled to lessen requirements
using Distributed
using LinearAlgebra
using Random

export contour_beyn

"""
    λv,V=contour_beyn([eltype,] nep;[tol,][displaylevel,][σ,][radius,][linsolvercreator,][quad_method,][N,][neigs,][k])

The function computes eigenvalues using Beyn's contour integral approach,
using an ellipse centered at `σ` with radii given in `radius`, or if only one `radius` is given,
the contour is a circle. The quadrature method is specified in `quad_method`
(`:ptrapz`,`ptrapz_parallel`, `:quadg`,`:quadg_parallel`,`:quadgk`). `k`
specifies the number of computed eigenvalues. `N` corresponds to the
number of quadrature points. Ellipses are the only supported contours. The
`linsolvercreator` must create a linsolver that can handle (rectangular) matrices
as right-hand sides, not only vectors. We integrate in complex arithmetic so
`eltype` must be complex type.

The kwargs `neigs` specifies the number of wanted eigvals, and `k` is the number
of columns in the matrix to be integrated (default `k=neigs+1`). If you give the `k` parameter and set `neigs=typemax(Int)` all found eigenvalues will be returned. The kwarg `sanity_check` decides if sorting and checking (and removal) of eigpairs should be done.
If disabled, the method
returns `k` (potentially inaccurate) eigpairs. The parameters `errmeasure` and
`tol` and `rank_drop_tol` are used for the sanity check,
to extract accurate eigenvalues.

# Example
```julia-repl
julia> using LinearAlgebra
julia> nep=nep_gallery("dep0");
julia> # Look for two eigvals in unit disk
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
                      radius::Union{Real,Tuple,Array}=1, # integration radius
                      quad_method::Symbol=:ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz
                      N::Integer=1000,  # Nof quadrature nodes
                      errmeasure::Function =
                      default_errmeasure(nep::NEP),
                      sanity_check=true,
                      rank_drop_tol=tol # Used in sanity checking
                      )where{T<:Number}

    # Geometry
    length(radius)==1 ? radius=(radius,radius) : nothing
    g(t) = complex(radius[1]*cos(t),radius[2]*sin(t)) # ellipse
    gp(t) = complex(-radius[1]*sin(t),radius[2]*cos(t)) # derivative

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
    Tv(λ) = local_linsolve(T(λ),Vh)
    f(t) = Tv(g(t))*gp(t)
    @ifd(print("Computing integrals"))


    local A0,A1
    if (quad_method == :quadg_parallel)
        @ifd(print(" using quadg_parallel"))
        (A0,A1)=quadg_parallel(f,g,0,2*pi,N);
    elseif (quad_method == :quadg)
        #@ifd(print(" using quadg"))
        #A0=quadg(f1,0,2*pi,N);
        #A1=quadg(f2,0,2*pi,N);
        error("disabled");
    elseif (quad_method == :ptrapz)
        @ifd(print(" using ptrapz"))
        (A0,A1)=ptrapz(f,g,0,2*pi,N);
    elseif (quad_method == :ptrapz_parallel)
        @ifd(print(" using ptrapz_parallel"))
        (A0,A1)=ptrapz_parallel(f,g,0,2*pi,N);
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
        sorted_index = sortperm(map(x->abs(σ-x), λ));
        return (λ[sorted_index],V[:,sorted_index])
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
    # and potentially remove eigenvalues if more than neigs.
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
       @warn "Rank-drop not detected, your eigvals may be correct, but the algorithm cannot verify. Try to increase k. This warning can be disabled with `sanity_check=false`." S
    end

    if (size(λgood,1)<neigs  && neigs < typemax(Int))
       @warn "We found less eigvals than requested. Try increasing domain, or decreasing `tol`. This warning can be disabled with `sanity_check=false`." S
    end

    return (λgood,Vgood)
end

#  Carries out Gauss quadrature (with N) discretization points
#  by call to @parallel
function quadg_parallel(f,g,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    S = @distributed (+) for i = 1:N
        temp = f(t[i])*w[i]
        [temp,temp*g(t[i])]
    end
    return S[1], S[2];
end

function quadg(f,g,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    S0 = zero(f(t[1])); S1 = zero(S0)
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    for i = 1:N
        temp = f(t[i])*w[i]
        S0 += temp
        S1 += temp*g(t[i])
    end
    return S0, S1;
end


# Trapezoidal rule for a periodic function f
function ptrapz(f,g,a,b,N)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S0 = zero(f(t[1])); S1 = zero(S0)
    for i = 1:N
        temp = f(t[i])
        S0 += temp
        S1 += temp*g(t[i])
    end
    return h*S0, h*S1;
end


function ptrapz_parallel(f,g,a,b,N)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S = @distributed (+) for i = 1:N
        temp = f(t[i])
        [temp,temp*g(t[i])]
    end
    return h*S[1], h*S[2];
end
