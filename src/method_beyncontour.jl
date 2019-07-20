# Beyn contour integral approach

#using QuadGK # Disabled to lessen requirements
using Distributed
using LinearAlgebra
using Random

export contour_beyn

"""
    λv,V=contour_beyn([eltype,] nep [,mintegrator];[tol,][logger,][σ,][radius,][linsolvercreator,][N,][neigs,][k])

The function computes eigenvalues using Beyn's contour integral approach,
using an ellipse centered at `σ` with radii given in `radius`, or if only one `radius` is given,
the contour is a circle. The numerical quadrature method is specified in `mintegrator`,
which is a type inheriting from `MatrixIntegrator`, by default
`MatrixTrapezoidal`. For a parallell implementation of the
integrator use `MatrixTrapezoidalParallel`.
 The integer `k`
specifies size of the probe subspace. `N` corresponds to the
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
julia> λv,V=contour_beyn(nep,radius=1,neigs=2);
julia> norm(compute_Mlincomb(nep,λv[1],V[:,1])) # Eigenpair 1
5.778617503485546e-15
julia> norm(compute_Mlincomb(nep,λv[2],V[:,2])) # Eigenpair 2
3.095638020248726e-14
```
# References
* Wolf-Jürgen Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications 436 (2012) 3839–3863
"""
contour_beyn(nep::NEP;params...)=contour_beyn(ComplexF64,nep;params...)
contour_beyn(nep::NEP,MIntegrator;params...)=contour_beyn(ComplexF64,nep,MIntegrator;params...)
function contour_beyn(::Type{T},
                      nep::NEP,
                      ::Type{MIntegrator}=MatrixTrapezoidal;
                      tol::Real=sqrt(eps(real(T))), # Note tol is quite high for this method
                      σ::Number=zero(complex(T)),
                      logger=0,
                      linsolvercreator::Function=backslash_linsolvercreator,
                      neigs::Integer=2, # Number of wanted eigvals
                      k::Integer=neigs+1, # Columns in matrix to integrate
                      radius::Union{Real,Tuple,Array}=1, # integration radius
                      N::Integer=1000,  # Nof quadrature nodes
                      errmeasure::ErrmeasureType = DefaultErrmeasure,
                      sanity_check=true,
                      rank_drop_tol=tol # Used in sanity checking
                      )where{T<:Number, MIntegrator<:MatrixIntegrator}

    @parse_logger_param!(logger)

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

    # Init errmeasure
    ermdata=init_errmeasure(errmeasure,nep);


    Random.seed!(10); # Reproducability
    Vh=Array{T,2}(randn(real(T),n,k)) # randn only works for real


    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end

    # Constructing integrands
    Tv(λ) = local_linsolve(T(λ),Vh)
    f(t) = Tv(g(t))*gp(t)
    push_info!(logger,"Computing integrals")


    local A0,A1
    AA=integrate_interval(MIntegrator, ComplexF64, f,
                        [ s-> one(s), g],0,2*pi,N,logger )
    A0=AA[:,:,1];
    A1=AA[:,:,2];


    # Don't forget scaling
    A0[:,:] = A0 ./(2im*pi);
    A1[:,:] = A1 ./(2im*pi);

    push_info!(logger,"Computing SVD prepare for eigenvalue extraction ",continues=true)
    V,S,W = svd(A0)
    p = count( S/S[1] .> rank_drop_tol);
    push_info!(logger," p=$p");

    V0 = V[:,1:p]
    W0 = W[:,1:p]
    B = (copy(V0')*A1*W0) * Diagonal(1 ./ S[1:p])

    # Extract eigenval and eigvec approximations according to
    # step 6 on page 3849 in the reference
    push_info!(logger,"Computing eigenvalues ")
    λ,VB=eigen(B)
    λ[:] = λ .+ σ
    push_info!(logger,"Computing eigenvectors ")
    V = V0 * VB;
    for i = 1:p
        normalize!(V[:,i]);
    end

    if (!sanity_check)
        sorted_index = sortperm(map(x->abs(σ-x), λ));
        inside_bool = (real(λ[sorted_index].-σ)/radius[1]).^2 + (imag(λ[sorted_index].-σ)/radius[2]).^2 .≤ 1
        inside_perm = sortperm(.!inside_bool)
        return (λ[sorted_index[inside_perm]],V[:,sorted_index[inside_perm]])
    end

    # Compute all the errors
    errmeasures=zeros(real(T),p);
    for i = 1:p
        errmeasures[i]=estimate_error(ermdata,λ[i],V[:,i]);
    end

    good_index=findall(errmeasures .< tol);

    # Index vector for sorted to distance from σ
    sorted_good_index=
       good_index[sortperm(map(x->abs(σ-x), λ[good_index]))];

    # check that eigenvalues are inside contour, move them to the end if they are outside
   inside_bool = (real(λ[sorted_good_index].-σ)/radius[1]).^2 + (imag(λ[sorted_good_index].-σ)/radius[2]).^2 .≤ 1
   if any(.!inside_bool)
       if neigs ≤ sum(inside_bool)
           @warn "found $(sum(.!inside_bool)) evals outside contour, $p inside. all $neigs returned evals inside, but possibly inaccurate. try increasing N, decreasing tol, or changing radius"
       else
           @warn "found $(sum(.!inside_bool)) evals outside contour, $p inside. last $(neigs-sum(inside_bool)) returned evals outside contour. possible inaccuracy. try increasing N, decreasing tol, or changing radius"
       end
   end
   sorted_good_inside_perm = sortperm(.!inside_bool)

    # Remove all eigpairs not sufficiently accurate
    # and potentially remove eigenvalues if more than neigs.
    local Vgood,λgood
    if( size(sorted_good_index,1) > neigs)
        found_evals=size(sorted_good_index,1);
        push_info!(logger,"Removing unwanted eigvals: neigs=$neigs<$found_evals=found_eigvals")
        Vgood=V[:,sorted_good_index[sorted_good_inside_perm][1:neigs]];
        λgood=λ[sorted_good_index[sorted_good_inside_perm][1:neigs]];
    else
        Vgood=V[:,sorted_good_index[sorted_good_inside_perm]];
        λgood=λ[sorted_good_index[sorted_good_inside_perm]];
    end

    if (p==k)
       @warn "Rank-drop not detected, your eigvals may be correct, but the algorithm cannot verify. Try to increase k. This warning can be disabled with `sanity_check=false`." S
    end

    if (size(λgood,1)<neigs  && neigs < typemax(Int))
       @warn "We found fewer eigvals than requested. Try increasing domain, or decreasing `tol`. This warning can be disabled with `sanity_check=false`." S
    end

    return (λgood,Vgood)
end
