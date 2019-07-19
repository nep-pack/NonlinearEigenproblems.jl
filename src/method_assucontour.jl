using Random;


export contour_block_SS

"""
    contour_block_SS

Higher order moments for contour integration (Asakura and Sakurai).


# Example

```julia-repl
julia> nep=SPMF_NEP([[0 1 ; 1 1.0], [1 0 ; 0 0]], [s->one(s),s->exp(1im*s^2)]);
julia> λ,V=contour_assu(nep,radius=3,neigs=6)
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
contour_block_SS(nep::NEP;params...)=contour_block_SS(ComplexF64,nep;params...);
function contour_block_SS(
    ::Type{T},
    nep::NEP;
    tol::Real=sqrt(eps(real(T))), # Note tol is quite high for this method
    σ::Number=zero(complex(T)),
    logger=0,
    linsolvercreator::Function=backslash_linsolvercreator,
    neigs::Integer=2, # Number of wanted eigvals (currently unused)
    k::Integer=2, # Columns in matrix to integrate
    radius::Union{Real,Tuple,Array}=1, # integration radius
    quad_method::Symbol=:ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz
    N::Integer=1000,  # Nof quadrature nodes
    K::Integer=3, # Nof moments
    L::Integer=3,
    errmeasure::ErrmeasureType = DefaultErrmeasure,
    sanity_check=true,
    rank_drop_tol=tol # Used in sanity checking
)where{T<:Number}

#    @parse_logger_param!(logger)


    n = size(nep,1);

    if (quad_method != :ptrapz)
        error("Only quad_method=:ptrapz currently supported")
    end

    Random.seed!(10); # Reproducability (not really)
    U = rand(T,n,L);
    V = rand(T,n,L);

    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
        print(".")
        local M0inv::LinSolver = linsolvercreator(nep,λ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end

    # The step-references refer to the JSIAM-paper

    # Quadrature points: So far only circle supported
    w = exp.(2im*pi*(0.5 .+ (0:(N-1)))/N);
    omega = σ .+ radius*w;

    # Step 2: Precompute all the linear systems
    FinvV =zeros(T,n,L,N);
    for k = 1:N
        FinvV[:,:,k]=local_linsolve(omega[k],V);
    end

    # Step 3-4: Compute all the integrals and store in Shat (

    Shat = zeros(T,n,L,2*K)
    Mhat = zeros(T,L,L,2*K)
    for k=0:(2*K-1)
        for j=0:N-1
            d=((omega[j+1]-σ)/radius)^(k+1)
            Shat[:,:,k+1] += d*FinvV[:,:,j+1]/N;
        end
        Mhat[:,:,k+1]=U'*Shat[:,:,k+1] # Step 4: Mhat=U'*Shat
    end


    # Construct H-matrices:
    m=K*L;
    Hhat=zeros(T,m,m)   # Hhat
    Hhat2=zeros(T,m,m)  # Hhat^{<}
    for i=1:K
        for j=1:K
            Hhat[(i-1)*L .+ (1:L), (j-1)*L .+ (1:L)] = Mhat[:,:,i+j-2+1];
            Hhat2[(i-1)*L .+ (1:L), (j-1)*L .+ (1:L)] = Mhat[:,:,i+j-1+1];
        end
    end


    # Extraction more similar to Algorithm 1
    # in https://arxiv.org/pdf/1510.02572.pdf

    F=svd(Hhat)
    UU=F.U;
    SS=F.S;
    VV=F.V

    # rank_drop_tol = δ in reference
    pp =   count( SS/SS[1] .> rank_drop_tol);
    mprime=pp; # To make closer to notation in paper



    # Pick relevant eigvecs
    UU_H1 = UU[:,1:mprime]
    VV_H1 = VV[:,1:mprime]


    # Step 7: Project the moment matrices:
    Hhat_mprime = UU_H1'*Hhat*VV_H1;
    Hhat2_mprime = UU_H1'*Hhat2*VV_H1;


    # Step 8:
    (xi,X)=eigen(Hhat2_mprime,Hhat_mprime)

    # Step 10: Extract eigpair

    # Compute S-matrix (by reshaping parts of Shat-tensor)
    S=zeros(T,n,L*K)
    for j=0:(K-1)
        S[:,j*L .+ (1:L)]=Shat[:,:,j+1]
    end
    V=S*VV_H1*X;

    λ=σ .+ radius*xi

    return λ,V
end
