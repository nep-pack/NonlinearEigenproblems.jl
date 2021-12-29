using Random;


export contour_block_SS
export just_test

"""
    contour_block_SS([eltype,] nep [,mintegrator];[tol,][logger,][σ,][radius,][linsolvercreator,][N,][neigs,][k,][L])

This is an implementation of the block_SS contour integral method which is
based on the computation of higher order moments.
The contour is an ellipse centered at `σ` with radii given in `radius`, or if only one `radius` is given, the contour is a circle. The numerical quadrature method is specified in `mintegrator`,
which is a type inheriting from `MatrixIntegrator`, by default
`MatrixTrapezoidal`. For a parallell implementation of the
integrator use `MatrixTrapezoidalParallel`.
 The integer `k`
specifies size of the probe subspace. `N` corresponds to the
number of quadrature points. The integer L specifies the number of moments.
Ellipses are the only supported contours. The
`linsolvercreator` must create a linsolver that can handle (rectangular) matrices
as right-hand sides, not only vectors. We integrate in complex arithmetic so
`eltype` must be complex type.



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
* Van Beeumen,  Meerbergen, Michiels. Connections between contour integration and rational Krylov methods for eigenvalue problems, 2016, TW673, https://lirias.kuleuven.be/retrieve/415487/
"""
contour_block_SS(nep::NEP;params...)=contour_block_SS(ComplexF64,nep;params...);
contour_block_SS(nep::NEP,MIntegrator;params...)=contour_block_SS(ComplexF64,nep,MIntegrator;params...);
function contour_block_SS(
    ::Type{T},
    nep::NEP,
    ::Type{MIntegrator}=MatrixTrapezoidal;
    tol::Real=sqrt(eps(real(T))), # Note tol is quite high for this method
    σ::Number=zero(complex(T)),
    logger=0,
    linsolvercreator=BackslashLinSolverCreator(),
    neigs=Inf, # Number of wanted eigvals (currently unused)
    k::Integer=3, # Columns in matrix to integrate
    radius::Union{Real,Tuple,Array}=1, # integration radius
    N::Integer=1000,  # Nof quadrature nodes
    K::Integer=3, # Nof moments
    errmeasure::ErrmeasureType = DefaultErrmeasure(nep),
    sanity_check=true,
    Shat_mode=:native, # native or JSIAM-mode
    rank_drop_tol=tol # Used in sanity checking
)where{T<:Number, MIntegrator<:MatrixIntegrator}

    @parse_logger_param!(logger)

    @printf "This is the HPC version"

    n = size(nep,1);


    # Notation: L in JSIAM-paper corresponds to k in Beyn's paper.
    # Input params the same as contourbeyn, but
    # the code is like JSIAM-paper
    L=k

    Random.seed!(10); # Reproducability (not really)
    U = rand(T,n,L);
    V = rand(T,n,L);

    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
        local M0inv::LinSolver = create_linsolver(linsolvercreator, nep, λ+σ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end

    # The step-references refer to the JSIAM-paper

    push_info!(logger,"Computing integrals")

    local Shat
    push_info!(logger,"Forming Mhat and Shat")
    Shat = zeros(T,n,L,2*K)
    Mhat = zeros(T,L,L,2*K)

    length(radius)==1 ? radius1=(radius,radius) : radius1=radius


    if (Shat_mode==:JSIAM)
        # This is the way the JSIAM-paper proposes to compute Shat
        if (length(radius)>1)
            error("JSIAM Shat_mode does not support ellipses");
        end

        # Quadrature points: Only circle supported
        w = exp.(2im*pi*(0.5 .+ (0:(N-1)))/N);
        omega = radius*w;

        push_info!(logger,"Forming all linear systems F(s)^{-1}V:",
                   continues=true)
        # Step 2: Precompute all the linear systems
        FinvV =zeros(T,n,L,N);
        for k = 1:N
            FinvV[:,:,k]=local_linsolve(omega[k],V);
        end
        push_info!(logger,"");

        # Step 3-4: Compute all the integrals and store in Shat (

        for k=0:(2*K-1)
            for j=0:N-1
                d=((omega[j+1])/radius)^(k+1)
                Shat[:,:,k+1] += d*FinvV[:,:,j+1]/N;
            end
        end

    elseif (Shat_mode==:native)

        # This deviates from the JSIAM-paper description, since
        # we do not precompute linear systems, but instead
        # compute linear system in combination with the quadrature.
        # It handles the scaling differently.
        # This version is also more extendable.

        g(t) = complex(radius1[1]*cos(t),radius1[2]*sin(t)) # ellipse
        gp(t) = complex(-radius1[1]*sin(t),radius1[2]*cos(t)) # derivative
        Tv(λ) = local_linsolve(T(λ),V)
        f(t) = Tv(g(t))*gp(t)/(2im*pi)

        gv=Vector{Function}(undef,2*K)
        for k=0:(2*K-1)
            gv[k+1]= s -> g(s)^k;
        end

        # Call the integrator
        Shat=integrate_interval(MIntegrator, ComplexF64,
                                f,gv,0,2*pi,N,logger )

    else
        error("Unknown Shat_mode: $Shat_mode")
    end

    for k=0:(2*K-1)
        Mhat[:,:,k+1]=U'*Shat[:,:,k+1] # Step 4: Mhat=U'*Shat
    end


    # Construct H-matrices:
    push_info!(logger,"Computing Hhat and Hhat^{<}")
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
    push_info!(logger,"Computing SVD prepare for eigenvalue extraction ",continues=true)
    F=svd(Hhat)
    UU=F.U;
    SS=F.S;
    VV=F.V

    # rank_drop_tol = δ in reference
    pp =   count( SS/SS[1] .> rank_drop_tol);
    mprime=pp; # To make closer to notation in paper

    push_info!(logger," mprime=$mprime");


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

    # Reverse the shift

    Shat_mode == :JSIAM  ? factor=radius : factor=1
    λ=σ .+ factor*xi

    return λ,V
end

function just_test()

    x = 10
    return x
end