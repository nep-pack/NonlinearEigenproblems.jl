using LinearAlgebra
using Random
using SparseArrays
using IterativeSolvers
using NonlinearEigenproblems

export nleigs_coefficients

"""
    nleigs(nep::NEP, Σ::AbstractVector{Complex{T}})

Find a few eigenvalues and eigenvectors of a nonlinear eigenvalue problem, using the `nleigs` algorithm.

# Arguments
- `nep`: An instance of a nonlinear eigenvalue problem.
- `Σ`: A vector containing the points of a polygonal target set in the complex plane.
- `Ξ`: A vector containing a discretization of the singularity set.
- `logger`: Level of display (0, 1, 2).
- `maxdgr`: Max degree of approximation.
- `minit`: Min number of iterations after linearization is converged.
- `maxit`: Max number of total iterations.
- `tollin`: Tolerance for convergence of linearization.
- `isfunm` : Whether to use matrix functions.
- `static`: Whether to use static version of NLEIGS.
- `leja`: Use of Leja-Bagby points (0 = no, 1 = only in expansion phase, 2 = always).
- `nodes`: Prefixed interpolation nodes (only when leja is 0 or 1).
- `blksize`: Block size for pre-allocation.

See [`augnewton`](@ref) for other parameters.


# Return values
- `λ`: Vector of eigenvalues of the nonlinear eigenvalue problem NLEP inside the target set Σ.
- `X`: Corresponding matrix of eigenvectors.
- `res`: Corresponding residuals.
- `details`: Solution details, if requested (see NleigsSolutionDetails).

# Example
```julia-repl
julia> nep=nep_gallery("dep0");
julia> unit_square = float([1+1im, 1-1im, -1-1im,-1+1im])
julia> (λ,v)=nleigs(nep,unit_square);
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
2.4522684986758914e-12
julia> norm(compute_Mlincomb(nep,λ[2],v[:,2]))
2.7572460495529512e-12
```
# References
- S. Guettel, R. Van Beeumen, K. Meerbergen, and W. Michiels. NLEIGS: A class
  of fully rational Krylov methods for nonlinear eigenvalue problems. SIAM J.
  Sci. Comput., 36(6), A2842-A2864, 2014.
- [NLEIGS Matlab toolbox](http://twr.cs.kuleuven.be/research/software/nleps/nleigs.php)
"""
nleigs_coefficients(nep, Σ; params...) = nleigs_coefficients(Float64, nep, Σ; params...)
function nleigs_coefficients(
        ::Type{T},
        nep::NEP,
        Σ::AbstractVector{CT}=Vector{CT}([-1.0-1im,-1+1im,+1+1im,1-1im]);
        Ξ::Vector{T} = [T(Inf)],
        logger = 0,
        maxdgr::Int = 100,
        maxit::Int = 200,
        tollin::T = 100*eps(T),
        isfunm::Bool = true,
        leja::Int = 1,
        nodes::Vector{CT} = Vector{CT}(),
        blksize::Int = 20,
        return_details::Bool = false
        ) where {T<:Real, CT<:Complex{T}}

#    @parse_logger_param!(logger)

    # The following variables are used when creating the return values, so put them in scope
    D = Vector{Matrix{CT}}()

    P = NonlinearEigenproblems.RKHelper.get_rk_nep(T, nep)
    n = size(nep, 1)
    n == 1 && (maxdgr = maxit + 1)
    computeD = (n <= 400) # for small problems, explicitly use generalized divided differences
    b = blksize


    nrmD = Array{T}(undef, 1)

    # Discretization of Σ --> Gamma & Leja-Bagby points
    if leja == 0 # use no leja nodes
        if isempty(nodes)
            error("Interpolation nodes must be provided via 'nodes' when no Leja-Bagby points ('leja' == 0) are used.")
        end
        gamma,_ = NonlinearEigenproblems.RKHelper.discretizepolygon(Σ)
        max_count = max(maxit,maxdgr)+2
        σ = repeat(reshape(nodes, :, 1), ceil(Int, max_count/length(nodes)), 1)
        _,ξ,β = NonlinearEigenproblems.RKHelper.lejabagby(σ[1:maxdgr+2], Ξ, gamma, maxdgr+2, true, P.p)
    elseif leja == 1 # use leja nodes in expansion phase
        if isempty(nodes)
            gamma,nodes = NonlinearEigenproblems.RKHelper.discretizepolygon(Σ, true)
        else
            gamma,_ = NonlinearEigenproblems.RKHelper.discretizepolygon(Σ)
        end
        nodes = repeat(reshape(nodes, :, 1), ceil(Int, (maxit+1)/length(nodes)), 1)
        σ,ξ,β = NonlinearEigenproblems.RKHelper.lejabagby(gamma, Ξ, gamma, maxdgr+2, false, P.p)
    else # use leja nodes in both phases
        gamma,_ = NonlinearEigenproblems.RKHelper.discretizepolygon(Σ)
        max_count =  max(maxit,maxdgr)+2
        σ,ξ,β = NonlinearEigenproblems.RKHelper.lejabagby(gamma, Ξ, gamma, max_count, false, P.p)
    end
    ξ[maxdgr+2] = NaN # not used
    if (!P.spmf || !isfunm) && length(σ) != length(unique(σ))
        error("All interpolation nodes must be distinct when no matrix " *
            "functions are used for computing the generalized divided differences.")
    end

    # Rational Newton coefficients
    range = 1:maxdgr+2
    if !P.spmf
        D = ratnewtoncoeffs(λ -> compute_Mder(nep, λ[1]), σ[range], ξ[range], β[range])
        nrmD[1] = norm(D[1]) # Frobenius norm
        sgdd = Matrix{CT}(undef, 0, 0)
    else
        # Compute scalar generalized divided differences
        sgdd = NonlinearEigenproblems.RKHelper.scgendivdiffs(σ[range], ξ[range], β[range], maxdgr, isfunm, get_fv(nep))
        # Construct first generalized divided difference
        computeD && push!(D, constructD(0, P, sgdd))
        # Norm of first generalized divided difference
        nrmD[1] = maximum(abs.(sgdd[:,1]))
    end
    if !isfinite(nrmD[1]) # check for NaN
        error("The generalized divided differences must be finite.");
    end

    # Rational Krylov
    expand = true
    kconv = trunc(Int, typemax(Int)/2)
    kn = n   # length of vectors in V
    l = 0    # number of vectors in V
    N = 0    # degree of approximations
    nbconv = 0 # number of converged lambdas inside Σ
    nblamin = 0 # number of lambdas inside Σ, converged or not
    kmax =  maxit
    k = 1
    while k <= kmax

        if expand

            # rational divided differences
            if P.spmf && computeD
                push!(D, constructD(k, P, sgdd))
            end
            N += 1

            # monitoring norms of divided difference matrices
            if !P.spmf
                push!(nrmD, norm(D[k+1])) # Frobenius norm
            else
                # The below can cause out of bounds in sgdd if there's
                # no convergence (also happens in MATLAB implementation)
                push!(nrmD, maximum(abs.(sgdd[:,k+1])))
            end
            if !isfinite(nrmD[k+1]) # check for NaN
                error("The generalized divided differences must be finite.");
            end
            if n > 1 && k >= 5 && k < kconv
                if sum(nrmD[k-3:k+1]) < 5*tollin
                    kconv = k - 1
                    expand = false
                    if leja == 1
                        # TODO: can we pre-allocate σ?
                        if length(σ) < kmax+1
                            resize!(σ, kmax+1)
                        end
                        σ[k+1:kmax+1] = nodes[1:kmax-k+1]
                    end
                    if !P.spmf || computeD
                        D = D[1:k]
                    end
                    ξ = ξ[1:k]
                    β = β[1:k]
                    nrmD = nrmD[1:k]

                    N -= 1
                    #NonlinearEigenproblems.push_info!(logger,
                    #           "Linearization converged after $kconv iterations")
                    #NonlinearEigenproblems.push_info!(logger,
                    #           "--> freeze linearization")
                elseif k == maxdgr+1
                    kconv = k
                    expand = false
                    if leja == 1
                        if length(σ) < kmax+1
                            resize!(σ, kmax+1)
                        end
                        σ[k+1:kmax+1] = nodes[1:kmax-k+1]
                    end

                    N -= 1
                    @warn "NLEIGS: Linearization not converged after $maxdgr iterations"
                    #NonlinearEigenproblems.push_info!(logger,
                    #           "--> freeze linearization")
                end
            end
        end

        l =  k


        # increment k
        k += 1
    end

    return D, β, ξ, σ
end

"Construct generalized divided difference for number `nb`."
function constructD(nb, P, sgdd::AbstractMatrix{CT}) where CT<:Complex{<:Real}
    n = size(P.nep, 1)
    BC = get_Av(P.nep)
    if !P.is_low_rank || nb <= P.p
        D = spzeros(CT, n, n)
        for ii = 1:(P.p+1+P.q)
            D += sgdd[ii,nb+1] * BC[ii]
        end
    else
        D = []
        for ii = 1:P.q
            d = sgdd[P.p+1+ii,nb+1] * P.L[ii]
            D = ii == 1 ? d : hcat(D, d)
        end
    end
    return D
end

"True for complex points `z` inside polygonal set `Σ`."
function in_Σ(z::AbstractVector{CT}, Σ::AbstractVector{CT}, tollin::T) where {T<:Real, CT<:Complex{T}}
    if length(Σ) == 2 && isreal(Σ)
        realΣ = real([Σ[1]; Σ[1]; Σ[2]; Σ[2]])
        imagΣ = [-tollin; tollin; tollin; -tollin]
    else
        realΣ = real(Σ)
        imagΣ = imag(Σ)
    end
    return map(p -> inpolygon(real(p), imag(p), realΣ, imagΣ), z)
end
