using NonlinearEigenproblems.RKHelper
using LinearAlgebra
using Random

export nleigs, NleigsSolutionDetails

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
- `tol`: Tolerance for residual.
- `tollin`: Tolerance for convergence of linearization.
- `v`: Starting vector.
- `errmeasure`: Function for error measure (residual norm). Called with arguments (λ,v).
- `isfunm` : Whether to use matrix functions.
- `static`: Whether to use static version of NLEIGS.
- `leja`: Use of Leja-Bagby points (0 = no, 1 = only in expansion phase, 2 = always).
- `nodes`: Prefixed interpolation nodes (only when leja is 0 or 1).
- `reusefact`: Reuse of matrix factorizations (0 = no, 1 = only after converged linearization, 2 = always).
- `blksize`: Block size for pre-allocation.
- `return_details`: Whether to return solution details (see NleigsSolutionDetails).
- `check_error_every`: Check for convergence / termination every this number of iterations.

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
5.028313023882492e-14
julia> norm(compute_Mlincomb(nep,λ[2],v[:,2]))
1.1937025845487509e-13
```
# References
- S. Guettel, R. Van Beeumen, K. Meerbergen, and W. Michiels. NLEIGS: A class
  of fully rational Krylov methods for nonlinear eigenvalue problems. SIAM J.
  Sci. Comput., 36(6), A2842-A2864, 2014.
- [NLEIGS Matlab toolbox](http://twr.cs.kuleuven.be/research/software/nleps/nleigs.php) (GPL License)
- [NLEIGS Matlab toolbox](https://bitbucket.org/roelvb/nleigs) (MIT License)
"""
nleigs(nep, Σ; params...) = nleigs(Float64, nep, Σ; params...)
function nleigs(
        ::Type{T},
        nep::NEP,
        Σ::AbstractVector{CT}=Vector{CT}([-1.0-1im,-1+1im,+1+1im,1-1im]);
        Ξ::Vector{T} = [T(Inf)],
        logger = 0,
        maxdgr::Int = 100,
        minit::Int = 20,
        maxit::Int = 200,
        linsolvercreator=DefaultLinSolverCreator(),
        tol::T = 1e-10,
        tollin::T = max(tol/10, 100*eps(T)),
        v::Vector{CT} = CT.(randn(T, size(nep, 1))),
        errmeasure::ErrmeasureType = ResidualErrmeasure(nep),
        isfunm::Bool = true,
        static::Bool = false,
        leja::Int = 1,
        nodes::Vector{CT} = Vector{CT}(),
        reusefact::Int = 1,
        blksize::Int = 20,
        return_details::Bool = false,
        check_error_every::Int = 5) where {T<:Real, CT<:Complex{T}}

    @parse_logger_param!(logger)

    # The following variables are used when creating the return values, so put them in scope
    D = Vector{Matrix{CT}}()
    conv = BitVector()
    lam = Vector{CT}()
    X = Matrix{CT}(undef, 0, 0)
    res = Vector{T}()

    P = get_rk_nep(T, nep)
    n = size(nep, 1)
    n == 1 && (maxdgr = maxit + 1)
    computeD = (n <= 400) # for small problems, explicitly use generalized divided differences
    b = blksize

    lin_solver_cache = LinSolverCache(CT, nep, linsolvercreator)

    # Initialization
    if static
        V = zeros(CT, n, 1)
    elseif !P.spmf
        V = zeros(CT, (b+1)*n, b+1)
    else
        if !P.is_low_rank || b < P.p
            V = zeros(CT, (b+1)*n, b+1)
        else
            V = zeros(CT, P.p*n+(b-P.p+1)*P.r, b+1)
        end
    end
    H = zeros(CT, b+1, b)
    K = zeros(CT, b+1, b)
    nrmD = Array{T}(undef, 1)
    if return_details
        Lam = zeros(CT, b, b)
        Res = zeros(T, b, b)
    end

    # Discretization of Σ --> Gamma & Leja-Bagby points
    if leja == 0 # use no leja nodes
        if isempty(nodes)
            error("Interpolation nodes must be provided via 'nodes' when no Leja-Bagby points ('leja' == 0) are used.")
        end
        gamma,_ = discretizepolygon(Σ)
        max_count = static ? maxit+maxdgr+2 : max(maxit,maxdgr)+2
        σ = repeat(reshape(nodes, :, 1), ceil(Int, max_count/length(nodes)), 1)
        _,ξ,β = lejabagby(σ[1:maxdgr+2], Ξ, gamma, maxdgr+2, true, P.p)
    elseif leja == 1 # use leja nodes in expansion phase
        if isempty(nodes)
            gamma,nodes = discretizepolygon(Σ, true)
        else
            gamma,_ = discretizepolygon(Σ)
        end
        nodes = repeat(reshape(nodes, :, 1), ceil(Int, (maxit+1)/length(nodes)), 1)
        σ,ξ,β = lejabagby(gamma, Ξ, gamma, maxdgr+2, false, P.p)
    else # use leja nodes in both phases
        gamma,_ = discretizepolygon(Σ)
        max_count = static ? maxit+maxdgr+2 : max(maxit,maxdgr)+2
        σ,ξ,β = lejabagby(gamma, Ξ, gamma, max_count, false, P.p)
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
        sgdd = scgendivdiffs(σ[range], ξ[range], β[range], maxdgr, isfunm, get_fv(nep))
        # Construct first generalized divided difference
        computeD && push!(D, constructD(0, P, sgdd))
        # Norm of first generalized divided difference
        nrmD[1] = maximum(abs.(sgdd[:,1]))
    end
    if !isfinite(nrmD[1]) # check for NaN
        error("The generalized divided differences must be finite.");
    end

    # Rational Krylov
    v = solve(lin_solver_cache, σ[1], v/norm(v), reusefact == 2)
    V[1:n,1] .= v ./ norm(v)
    expand = true
    kconv = trunc(Int, typemax(Int)/2)
    kn = n   # length of vectors in V
    l = 0    # number of vectors in V
    N = 0    # degree of approximations
    nbconv = 0 # number of converged lambdas inside Σ
    nblamin = 0 # number of lambdas inside Σ, converged or not
    kmax = static ? maxit + maxdgr : maxit
    k = 1
    while k <= kmax
        # resize matrices if we're starting a new block
        if l > 0 && (b == 1 || mod(l+1, b) == 1)
            nb = round(Int, 1 + l/b)
            Vrows = size(V, 1)
            if !P.spmf
                Vrows = kn+b*n
            else
                if expand
                    if !P.is_low_rank || l + b < P.p
                        Vrows = kn+b*n
                    elseif l < P.p-1
                        Vrows = P.p*n+(nb*b-P.p+1)*P.r
                    else # l >= P.p-1
                        Vrows = kn+b*P.r
                    end
                end
            end
            V = resize_matrix(V, Vrows, nb*b+1)
            H = resize_matrix(H, size(H, 1) + b, size(H, 2) + b)
            K = resize_matrix(K, size(K, 1) + b, size(K, 2) + b)
            if return_details
                Lam = resize_matrix(Lam, size(Lam, 1) + b, size(Lam, 2) + b)
                Res = resize_matrix(Res, size(Res, 1) + b, size(Res, 2) + b)
            end
        end

        if expand
            # set length of vectors
            if !P.is_low_rank || k < P.p
                kn += n
            else
                kn += P.r
            end

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
                    if static
                        kmax = maxit + kconv
                    end
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
                    if static
                        if !P.is_low_rank || k < P.p
                            kn -= n
                        else
                            kn -= P.r
                        end
                        V = resize_matrix(V, kn, b+1)
                    end
                    N -= 1
                    push_info!(logger,
                               "Linearization converged after $kconv iterations")
                    push_info!(logger,
                               "--> freeze linearization")
                elseif k == maxdgr+1
                    kconv = k
                    expand = false
                    if leja == 1
                        if length(σ) < kmax+1
                            resize!(σ, kmax+1)
                        end
                        σ[k+1:kmax+1] = nodes[1:kmax-k+1]
                    end
                    if static
                        V = resize_matrix(V, kn, b+1)
                    end
                    N -= 1
                    @warn "NLEIGS: Linearization not converged after $maxdgr iterations"
                    push_info!(logger,
                               "--> freeze linearization")
                end
            end
        end

        l = static ? k - N : k

        if !static || (static && !expand)
            # shift-and-invert
            t = [zeros(l-1); 1]    # continuation combination
            wc = V[1:kn, l]        # continuation vector
            w = backslash(wc, P, lin_solver_cache, reusefact, computeD, σ, k, D, β, N, ξ, expand, kconv, sgdd)

            # orthogonalization
            Vview = view(V, 1:kn, 1:l)
            H[l+1,l] = orthogonalize_and_normalize!(Vview, w, view(H, 1:l,l), DGKS())
            K[1:l,l] .= view(H, 1:l, l) .* σ[k+1] .+ t
            K[l+1,l] = H[l+1,l] * σ[k+1]
            V[1:kn,l+1] = w
#            @printf("new vector V: size = %s, sum = %s\n", size(V[1:kn,l+1]), sum(sum(V[1:kn,l+1])))
        end

        function check_convergence(all)
            lambda, S = eigen(K[1:l,1:l], H[1:l,1:l])

            # select eigenvalues
            if !all
                lamin = in_Σ(lambda, Σ, tol)
                ilam = [1:l;][lamin]
                lam = lambda[ilam]
            else
                ilam = [1:l;][isfinite.(lambda)]
                lam = lambda[ilam]
                lamin = in_Σ(lam, Σ, tol)
            end

            nblamin = sum(lamin)
            for i = ilam
                S[:,i] /= norm(H[1:l+1,1:l] * S[:,i])
            end
            X = V[1:n,1:l+1] * (H[1:l+1,1:l] * S[:,ilam])
            for i = 1:size(X,2)
                X[:,i] /= norm(X[:,i])
            end

            # compute residuals & check for convergence
            res = map(i -> estimate_error(errmeasure,lam[i], X[:,i]), 1:length(lam))
            conv = abs.(res) .< tol
            if all
                resall = fill(T(NaN), l, 1)
                resall[ilam] = res
                # sort complex numbers by magnitude, then angle
                si = sortperm(lambda, lt = (a,b) -> abs(a) < abs(b) || (abs(a) == abs(b) && angle(a) < angle(b)))
                Res[1:l,l] = resall[si]
                Lam[1:l,l] = lambda[si]
                conv .&= lamin
            end

            nbconv = isempty(conv) ? 0 : sum(conv)

            iteration = static ? k - N : k
            push_info!(logger,
                       "  iteration $iteration: $nbconv of $nblamin < $tol")
        end

        # Ritz pairs
        if !return_details && (
            (!expand && k >= N + minit && mod(k-(N+minit), check_error_every) == 0) ||
            (k >= kconv + minit && mod(k-(kconv+minit), check_error_every) == 0) || k == kmax)
            check_convergence(false)
        elseif return_details && (!static || (static && !expand))
            check_convergence(true)
        end

        # stopping
        if ((!expand && k >= N + minit) || k >= kconv + minit) && nblamin == nbconv
            break
        end

        # increment k
        k += 1
    end

    details = NleigsSolutionDetails{T,CT}()

    if return_details
        Lam = Lam[1:l,1:l]
        Res = Res[1:l,1:l]
        σ = σ[1:k]
        if expand
            ξ = ξ[1:k]
            β = β[1:k]
            nrmD = nrmD[1:k]
            @warn "NLEIGS: Linearization not converged after $maxdgr iterations"
        end
        details = NleigsSolutionDetails(Lam, Res, σ, ξ, β, nrmD, kconv)
    end

    return lam[conv], X[:,conv], res[conv], details
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

"Backslash or left matrix divide for continuation vector `wc`."
function backslash(wc, P, lin_solver_cache, reusefact, computeD, σ, k, D, β, N, ξ, expand, kconv, sgdd)
    n = size(P.nep, 1)::Int
    shift = σ[k+1]

    # construction of B*wc
    Bw = zeros(eltype(wc), size(wc))
    # first block (if low rank)
    if P.is_low_rank
        i0b = (P.p-1)*n + 1
        i0e = P.p*n
        if !P.spmf || computeD
            Bw[1:n] = -D[P.p+1] * wc[i0b:i0e] / β[P.p+1]
        else
            Bw[1:n] = -sum(reshape(P.BBCC * wc[i0b:i0e], n, :) .* transpose(sgdd[:,P.p+1]), dims = 2) / β[P.p+1];
        end
    end
    # other blocks
    i0b = 1
    i0e = n
    for ii = 1:N
        # range of block i+1
        i1b = i0e + 1
        if !P.is_low_rank || ii < P.p
            i1e = i0e + n
        else
            i1e = i0e + P.r
        end
        # compute block i+1
        if !P.is_low_rank || ii != P.p
            Bw[i1b:i1e] = wc[i0b:i0e] + β[ii+1]/ξ[ii]*wc[i1b:i1e]
        else
            Bw[i1b:i1e] = P.UU' * wc[i0b:i0e] + β[ii+1]/ξ[ii]*wc[i1b:i1e]
        end
        # range of next block i
        i0b = i1b
        i0e = i1e
    end

    # construction of z0
    z = copy(Bw)
    i1b = n + 1
    if !P.is_low_rank || P.p > 1
        i1e = 2*n
    else
        i1e = n + P.r
    end
    nu = β[2]*(1 - shift/ξ[1])
    z[i1b:i1e] = 1/nu * z[i1b:i1e]
    for ii = 1:N
        # range of block i+2
        i2b = i1e + 1
        if !P.is_low_rank || ii < P.p-1
            i2e = i1e + n
        else
            i2e = i1e + P.r
        end
        # add extra term to z0
        if !P.spmf || computeD
            if !P.is_low_rank || ii != P.p
                z[1:n] -= D[ii+1] * z[i1b:i1e]
            end
        else
            if !P.is_low_rank || ii < P.p
                z[1:n] -= sum(reshape(P.BBCC * z[i1b:i1e], n, :) .* transpose(sgdd[:,ii+1]), dims = 2)
            elseif ii > P.p
                dd = sgdd[P.p+2:end,ii+1]
                @inbounds for i = 1:length(P.iLr)
                    for j = 1:length(P.LL[i].nzval)
                        ind = P.LL[i].nzind[j]
                        z[P.iLr[i]] -= P.LL[i].nzval[j] * z[i1b+ind-1] * dd[P.iL[ind]]
                    end
                end
            end
        end
        # update block i+2
        if ii < N
            mu = shift - σ[ii+1]
            nu = β[ii+2] * (1 - shift/ξ[ii+1])
            if !P.is_low_rank || ii != P.p-1
                z[i2b:i2e] = 1/nu * z[i2b:i2e] + mu/nu * z[i1b:i1e]
            else # i == p-1
                z[i2b:i2e] = 1/nu * z[i2b:i2e] + mu/nu * P.UU'*z[i1b:i1e]
            end
        end
        # range of next block i+1
        i1b = i2b
        i1e = i2e
    end

    # solving Alam x0 = z0
    w = zeros(eltype(wc), size(wc))
    add_to_cache = ((!expand || k > kconv) && reusefact == 1) || reusefact == 2
    w[1:n] = solve(lin_solver_cache, shift, z[1:n]/β[1], add_to_cache)

    # substitutions x[i+1] = mu/nu*x[i] + 1/nu*Bw[i+1]
    i0b = 1
    i0e = n
    for ii = 1:N
        # range of block i+1
        i1b = i0e + 1
        if !P.is_low_rank || ii < P.p
            i1e = i0e + n
        else
            i1e = i0e + P.r
        end
        # compute block i+1
        mu = shift - σ[ii]
        nu = β[ii+1] * (1 - shift/ξ[ii])
        if !P.is_low_rank || ii != P.p
            w[i1b:i1e] = mu/nu * w[i0b:i0e] + 1/nu * Bw[i1b:i1e]
        else
            w[i1b:i1e] = mu/nu * P.UU'*w[i0b:i0e] + 1/nu * Bw[i1b:i1e]
        end
        # range of next block i
        i0b = i1b
        i0e = i1e
    end

    return w
end

"True for complex points `z` inside polygonal set `Σ`."
function in_Σ(z::AbstractVector{CT}, Σ::AbstractVector{CT}, tol::T) where {T<:Real, CT<:Complex{T}}
    if length(Σ) == 2 && isreal(Σ)
        realΣ = real([Σ[1]; Σ[1]; Σ[2]; Σ[2]])
        imagΣ = [-tol; tol; tol; -tol]
    else
        realΣ = real(Σ)
        imagΣ = imag(Σ)
    end
    return map(p -> inpolygon(real(p), imag(p), realΣ, imagΣ), z)
end

function resize_matrix(A, rows, cols)
    resized = zeros(eltype(A), rows, cols)
    resized[1:size(A, 1), 1:size(A, 2)] = A
    return resized
end

struct NleigsSolutionDetails{T<:Real, CT<:Complex{T}}
    "matrix of Ritz values in each iteration"
    Lam::AbstractMatrix{CT}

    "matrix of residuals in each iteraion"
    Res::AbstractMatrix{T}

    "vector of interpolation nodes"
    σ::AbstractVector{CT}

    "vector of poles"
    ξ::AbstractVector{T}

    "vector of scaling parameters"
    β::AbstractVector{T}

    "vector of norms of generalized divided differences (in function handle
    case) or maximum of absolute values of scalar divided differences in
    each iteration (in matrix function case)"
    nrmD::AbstractVector{T}

    "number of iterations until linearization converged"
    kconv::Int
end

NleigsSolutionDetails{T,CT}() where {T<:Real, CT<:Complex{T}} = NleigsSolutionDetails(
    Matrix{CT}(undef, 0, 0), Matrix{T}(undef, 0, 0), Vector{CT}(),
    Vector{T}(), Vector{T}(), Vector{T}(), 0)
