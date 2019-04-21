using Random;

export contour_assu

"""
    contour_assu

Higher order moments for contour integration (Asakura and Sakurai).


# Example

```julia-repl
julia> nep=SPMF_NEP([[0 1 ; 1 1.0], [1 0 ; 0 0]], [s->one(s),s->exp(1im*s^2)]);
julia> λ,V=countour_assu(nep,radius=3,neigs=6)
julia> @show λ
6-element Array{Complex{Float64},1}:
  4.496403249731884e-15 + 2.506628274630998im
        -2.506628274631 - 2.8727020762175925e-15im
  3.219972424519104e-16 - 2.5066282746310034im
     2.5066282746310096 - 1.1438072192922029e-15im
 -2.3814273710772784e-7 - 7.748469160458366e-8im
   2.381427350935646e-7 + 7.748467479992284e-8im
```

# References

* Asakura, Sakurai, Tadano, Ikegami, Kimura, A numerical method for nonlinear eigenvalue problems using contour integrals, JSIAM Letters, 2009 Volume 1 Pages 52-55
* Algorithm 2: Wolf-Jürgen Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications 436 (2012) 3839–3863
* Figure 5.3: Tisseur, Güttel. Nonlinear eigenvalue problems, 2017, vol 26
* Van Beeumen,  Meerbergen, Michiels. Connections between contour integration and rational Krylov methods for eigenvalue problems, 2016, TW673, https://lirias.kuleuven.be/retrieve/415487/
"""
contour_assu(nep::NEP;params...)=contour_assu(ComplexF64,nep;params...)
function countour_assu(
    ::Type{T},
    nep::NEP;
    tol::Real=sqrt(eps(real(T))), # Note tol is quite high for this method
    σ::Number=zero(complex(T)),
    displaylevel::Integer=0,
    linsolvercreator::Function=backslash_linsolvercreator,
    neigs::Integer=2, # Number of wanted eigvals
    k::Integer=2, # Columns in matrix to integrate
    radius::Union{Real,Tuple,Array}=1, # integration radius
    quad_method::Symbol=:ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz
    N::Integer=1000,  # Nof quadrature nodes
    errmeasure::ErrmeasureType = DefaultErrmeasure,
    sanity_check=true,
    rank_drop_tol=tol # Used in sanity checking
)where{T<:Number}


    n = size(nep,1);

    # We restrict to the block-SS version so ell = r.
    ell = k;
    r = k;

    Random.seed!(10); # Reproducability
    L = rand(T,n,ell);
    R = rand(T,n,r);


    @show radius
    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
        @ifd(print("."))
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end


    # Quadrature points: So far only circle supported
    w = exp.(2im*pi*(1:N)/N);
    z = σ .+ radius*w;

    # Precompute all the integrals (last index is p)
    # (Use a tensor for precomputation as in Tisseur &  Guettel Figure 5.3)
    A =zeros(T,ell,r,2*neigs);
    for k = 1:N
        Fz=L'*local_linsolve(z[k],R);
        for j = 0:2*neigs-1
            A[:,:,j+1] += (w[k]^j*radius*w[k]/N)*Fz;
        end
    end

    @show size(A)
    @show neigs
    # Compute block matrices: B0 and B1:
    B0 = zeros(T,neigs*ell,neigs*r); B1 = copy(B0);
    for row = 0:neigs-1
        for col = 0:neigs-1
            B0[row*ell  .+ (1:ell),  col*r .+ (1:r)]=A[:,:,row+col+1]
            B1[row*ell  .+ (1:ell),  col*r .+ (1:r)]=A[:,:,row+col+2]
        end
    end
    # Use B0 for eigenvalue extraction
    (V,S,W) = svd(B0);

    @show neigs
    @show size(V)
    # Extract active subspace

    p = count( S/S[1] .> rank_drop_tol);

    V0 = V[:,1:p];
    W0 = W[:,1:p]
    BB = (copy(V0')*B1*W0) * Diagonal(1 ./ S[1:p])
    λ,VB = eigen(BB);
    λ[:] = radius*λ .+ σ


    # Still not sure how to extract eigevectors in the best way
    V=VB # Incorrect


    return λ,V
end
