include("lusolver.jl")
include("discretizepolygon.jl")
include("lejabagby.jl")
include("ratnewtoncoeffs.jl")
include("ratnewtoncoeffsm.jl")

#=
NLEIGS  Find a few eigenvalues and eigenvectors of a NLEP
   lambda = NLEIGS(NLEP,Sigma,[Xi]) returns a vector of eigenvalues of the
   nonlinear eigenvalue problem NLEP inside the target set Sigma. NLEP is a
   structure representing the nonlinear eigenvalue problem as a function
   handle:
     NLEP.Fun:  function handle A = NLEP.Fun(lambda)
     NLEP.n:    size of A
   or as a sum of constant matrices times scalar functions:
     NLEP.B:    array of matrices of the polynomial part
     NLEP.C:    array of matrices of the nonlinear part
     NLEP.L:    array of L-factors of matrices C (optional)
     NLEP.U:    array of U-factors of matrices C (optional)
     NLEP.f:    array of the nonlinear scalar or matrix functions
   Sigma is a vector containing the points of a polygonal target set in the
   complex plane. Xi is a vector containing a discretization of the singularity
   set. In case the singularity set is omitted Xi is set to Inf.

   [X,lambda,res] = NLEIGS(NLEP,Sigma,[Xi]) returns a matrix X of eigenvectors
   and a vector of corresponding eigenvalues of the nonlinear eigenvalue
   problem NLEP inside the target set Sigma. The vector res contains the
   corresponding residuals.

   [X,lambda,res] = NLEIGS(NLEP,Sigma,Xi,options) sets the algorithm's
   parameters to the values in the structure options:
     options.disp:     level of display
                       [ {0} | 1 | 2 ]
     options.maxdgr:   max degree of approximation
                       [ positive integer {100} ]
     options.minit:    min number of iter. after linearization is converged
                       [ positive integer {20} ]
     options.maxit:    max number of total iterations
                       [ positive integer {200} ]
     options.tolres:   tolerance for residual
                       [ positive scalar {1e-10} ]
     options.tollin:   tolerance for convergence of linearization
                       [ positive scalar {1e-11} ]
     options.v0:       starting vector
                       [ vector array {randn(n,1)} ]
     options.funres:   function handle for residual (Lambda[vector], X[matrix])
                       [ @(Lambda,X) {norm(A(lam)*x)} ]
     options.isfunm:   use matrix functions
                       [ boolean {true} ]
     options.static:   static version of nleigs
                       [ boolean {false} ]
     options.leja:     use of Leja-Bagby points
                       [ 0 (no) | {1} (only in expansion phase) | 2 (always) ]
     options.nodes:    prefixed interpolation nodes (only with leja is 0 or 1)
                       [ vector array | {[]} ]
     options.reuselu:  reuse of LU-factorizations of A(sigma)
                       [ 0 (no) | {1} (only after conv. lin) | 2 (always) ]
     options.blksize:  block size for pre-allocation
                       [ positive integer {20} ]

   [X,lambda,res,info] = NLEIGS(NLEP,Sigma,[Xi]) also returns a structure info:
     info.Lam:    matrix of Ritz values in each iteration
     info.Res:    matrix of residuals in each iteraion
     info.sigma:  vector of interpolation nodes
     info.xi:     vector of poles
     info.beta:   vector of scaling parameters
     info.nrmD:   vector of norms of generalized divided differences (in
                  function handle case) or maximum of absolute values of
                  scalar divided differences in each iteration (in matrix
                  function case)

   See also EIG, EIGS, NLEIGSPLOT.

   Reference:
   S. Guettel, R. Van Beeumen, K. Meerbergen, and W. Michiels. NLEIGS: A class
   of fully rational Krylov methods for nonlinear eigenvalue problems. SIAM J.
   Sci. Comput., 36(6), A2842âA2864, 2014.

   Roel Van Beeumen
   April 5, 2016
=#
function nleigs(nep::NEP, Sigma::AbstractVector{Complex{T}}; Xi::AbstractVector{T} = [Inf], options::Dict = Dict(), return_info = false) where T<:Real
    # The following variables are used when creating the return values, so put them in scope
    D = Vector{Matrix{Complex{T}}}(0)
    conv = BitVector(0)
    lam = Vector{Complex{T}}(0)
    X = Matrix{Complex{T}}(0, 0)
    res = Vector{T}(0)

    element_type = eltype(Sigma)

    #@code_warntype prepare_inputs(nep, Sigma, Xi, options)
    #return
    spmf,iL,L,LL,UU,p,q,r,Sigma,leja,nodes,Xi,tollin,
        tolres,maxdgr,minit,maxit,isfunm,static,v0,reuselu,funres,b,computeD,
        resfreq,verbose,BBCC = prepare_inputs(nep, Sigma, Xi, options)
    n = nep.n

    # Initialization
    if static
        V = zeros(element_type, n, 1)
    elseif !spmf
        V = zeros(element_type, (b+1)*n, b+1)
    else
        if r == 0 || b < p
            V = zeros(element_type, (b+1)*n, b+1)
        else
            V = zeros(element_type, p*n+(b-p+1)*r, b+1)
        end
    end
    H = zeros(element_type, b+1, b)
    K = zeros(element_type, b+1, b)
    nrmD = Array{Float64}(1)
    if return_info
        Lam = zeros(element_type, b, b)
        Res = zeros(b, b)
    end

    # Discretization of Sigma --> Gamma & Leja-Bagby points
    if leja == 0 # use no leja nodes
        if isempty(nodes)
            error("Interpolation nodes must be provided via 'options[\"nodes\"]' when no Leja-Bagby points (options[\"leja\"] == 0) are used.")
        end
        gamma,_ = discretizepolygon(Sigma)
        max_count = static ? maxit+maxdgr+2 : max(maxit,maxdgr)+2
        sigma = repmat(reshape(nodes, :, 1), ceil(Int, max_count/length(nodes)), 1)
        _,xi,beta = lejabagby(sigma[1:maxdgr+2], Xi, gamma, maxdgr+2, true, p)
    elseif leja == 1 # use leja nodes in expansion phase
        if isempty(nodes)
            gamma,nodes = discretizepolygon(Sigma, true)
        else
            gamma,_ = discretizepolygon(Sigma)
        end
        nodes = repmat(reshape(nodes, :, 1), ceil(Int, (maxit+1)/length(nodes)), 1)
        sigma,xi,beta = lejabagby(gamma, Xi, gamma, maxdgr+2, false, p)
    else # use leja nodes in both phases
        gamma,_ = discretizepolygon(Sigma)
        max_count = static ? maxit+maxdgr+2 : max(maxit,maxdgr)+2
        sigma,xi,beta = lejabagby(gamma, Xi, gamma, max_count, false, p)
    end
    xi[maxdgr+2] = NaN # not used
    if (!spmf || !isfunm) && length(sigma) != length(unique(sigma))
        error("All interpolation nodes must be distinct when no matrix " *
            "functions are used for computing the generalized divided differences.")
    end

    # Rational Newton coefficients
    range = 1:maxdgr+2
    if !spmf
        D = ratnewtoncoeffs(λ -> compute_Mder(nep, λ[1]), sigma[range], xi[range], beta[range])
        nrmD[1] = vecnorm(D[1]) # Frobenius norm
        sgdd = Matrix{Complex{T}}(0, 0)
    else
        # Compute scalar generalized divided differences
        sgdd = scgendivdiffs(sigma[range], xi[range], beta[range], maxdgr, isfunm, nep.spmf.fi)
        # Construct first generalized divided difference
        computeD && push!(D, constructD(0, L, n, p, q, r, nep.spmf.A, sgdd))
        # Norm of first generalized divided difference
        nrmD[1] = maximum(abs.(sgdd[:,1]))
    end
    if !isfinite(nrmD[1]) # check for NaN
        error("The generalized divided differences must be finite.");
    end
    lureset()

    # Rational Krylov
    if reuselu == 2
        v0 = lusolve(λ -> compute_Mder(nep, λ), sigma[1], v0/norm(v0))
    else
        v0 = compute_Mder(nep, sigma[1]) \ (v0/norm(v0))
    end
    V[1:n,1] .= v0 ./ norm(v0)
    expand = true
    kconv = Inf
    kn = n   # length of vectors in V
    l = 0    # number of vectors in V
    N = 0    # degree of approximations
    nbconv = 0 # number of converged lambdas inside sigma
    nblamin = 0 # number of lambdas inside sigma, converged or not
    kmax = static ? maxit + maxdgr : maxit
    k = 1
    while k <= kmax
        # resize matrices if we're starting a new block
        if l > 0 && (b == 1 || mod(l+1, b) == 1)
            nb = round(Int, 1 + l/b)
            Vrows = size(V, 1)
            if !spmf
                Vrows = kn+b*n
            else
                if expand
                    if r == 0 || l + b < p
                        Vrows = kn+b*n
                    elseif l < p-1
                        Vrows = p*n+(nb*b-p+1)*r
                    else # l >= p-1
                        Vrows = kn+b*r
                    end
                end
            end
            V = resize_matrix(V, Vrows, nb*b+1)
            H = resize_matrix(H, size(H, 1) + b, size(H, 2) + b)
            K = resize_matrix(K, size(K, 1) + b, size(K, 2) + b)
            if return_info
                Lam = resize_matrix(Lam, size(Lam, 1) + b, size(Lam, 2) + b)
                Res = resize_matrix(Res, size(Res, 1) + b, size(Res, 2) + b)
            end
        end

        if expand
            # set length of vectors
            if r == 0 || k < p
                kn += n
            else
                kn += r
            end

            # rational divided differences
            if spmf && computeD
                push!(D, constructD(k, L, n, p, q, r, nep.spmf.A, sgdd))
            end
            N += 1

            # monitoring norms of divided difference matrices
            if !spmf
                push!(nrmD, vecnorm(D[k+1])) # Frobenius norm
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
                        # TODO: can we pre-allocate sigma?
                        if length(sigma) < kmax+1
                            resize!(sigma, kmax+1)
                        end
                        sigma[k+1:kmax+1] = nodes[1:kmax-k+1]
                    end
                    if !spmf || computeD
                        D = D[1:k]
                    end
                    xi = xi[1:k]
                    beta = beta[1:k]
                    nrmD = nrmD[1:k]
                    if static
                        if r == 0 || k < p
                            kn -= n
                        else
                            kn -= r
                        end
                        V = resize_matrix(V, kn, b+1)
                    end
                    N -= 1
                    if verbose > 0
                        println("Linearization converged after $kconv iterations")
                        println(" --> freeze linearization")
                    end
                elseif k == maxdgr+1
                    kconv = k
                    expand = false
                    if leja == 1
                        if length(sigma) < kmax+1
                            resize!(sigma, kmax+1)
                        end
                        sigma[k+1:kmax+1] = nodes[1:kmax-k+1]
                    end
                    if static
                        V = resize_matrix(V, kn, b+1)
                    end
                    N -= 1
                    warn("NLEIGS:noconvergence: Linearization not converged after $maxdgr iterations")
                    if verbose > 0
                        println(" --> freeze linearization")
                    end
                end
            end
        end

        l = static ? k - N : k

        if !static || (static && !expand)
            # shift-and-invert
            t = [zeros(l-1); 1]    # continuation combination
            wc = V[1:kn, l]        # continuation vector
            w = backslash(wc, nep, spmf, iL, LL, UU, p, q, r, reuselu, computeD, sigma, k, D, beta, N, xi, expand, kconv, BBCC, sgdd)

            # orthogonalization
            normw = norm(w)
            Vview = view(V, 1:kn, 1:l)
            h = Vview' * w
            w .-= Vview * h
            H[1:l,l] = h
            eta = 1/sqrt(2)       # reorthogonalization constant
            if norm(w) < eta * normw
                h = Vview' * w
                w .-= Vview * h
                H[1:l,l] .+= h
            end
            K[1:l,l] .= view(H, 1:l, l) .* sigma[k+1] .+ t

            # new vector
            H[l+1,l] = norm(w)
            K[l+1,l] = H[l+1,l] * sigma[k+1]
            V[1:kn,l+1] .= w ./ H[l+1,l]
#            @printf("new vector V: size = %s, sum = %s\n", size(V[1:kn,l+1]), sum(sum(V[1:kn,l+1])))
        end

        function update_lambdas(all)
            lambda, S = eig(K[1:l,1:l], H[1:l,1:l])

            # select eigenvalues
            if !all
                lamin = in_sigma(lambda, Sigma, tolres)
                ilam = [1:l;][lamin]
                lam = lambda[ilam]
            else
                ilam = [1:l;][isfinite.(lambda)]
                lam = lambda[ilam]
                lamin = in_sigma(lam, Sigma, tolres)
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
            res = funres(lam, X) # needed for residual: funA, Ahandle, BB, pf, CC, f, p, q
            conv = abs.(res) .< tolres
            if all
                resall = fill(NaN, l, 1)
                resall[ilam] = res
                # sort complex numbers by magnitude, then angle
                si = sortperm(lambda, lt = (a,b) -> abs(a) < abs(b) || (abs(a) == abs(b) && angle(a) < angle(b)))
                Res[1:l,l] = resall[si]
                Lam[1:l,l] = lambda[si]
                conv .&= lamin
            end

            nbconv = sum(conv)
            if verbose > 0
                iteration = static ? k - N : k
                println("  iteration $iteration: $nbconv of $nblamin < $tolres")
            end
        end

        # Ritz pairs
        if !return_info && (
            (!expand && k >= N + minit && mod(k-(N+minit), resfreq) == 0) ||
            (k >= kconv + minit && mod(k-(kconv+minit), resfreq) == 0) || k == kmax)
            update_lambdas(false)
        elseif return_info && (!static || (static && !expand))
            update_lambdas(true)
        end

        # stopping
        if ((!expand && k >= N + minit) || k >= kconv + minit) && nblamin == nbconv
            break
        end

        # increment k
        k += 1
    end
    lureset()

    info = Dict()

    if return_info
        Lam = Lam[1:l,1:l]
        Res = Res[1:l,1:l]
        sigma = sigma[1:k]
        if expand
            xi = xi[1:k]
            beta = beta[1:k]
            nrmD = nrmD[1:k]
            warn("NLEIGS:noconvergence: Linearization not converged after $maxdgr iterations")
        end
        info = Dict(
            "Lam" => Lam,
            "Res" => Res,
            "sigma" => sigma,
            "xi" => xi,
            "beta" => beta,
            "nrmD" => nrmD,
            "kconv" => kconv)
    end

    return X[:,conv], lam[conv], res[conv], info
end

# checkInputs: error checks the inputs to NLEP and also derives some variables from them:
#
#   funA       function handle for A(lambda)
#   Ahandle    is true if NLEP is given as function handle
#   B          cell array {B0,B1,...,Bp}
#   BB         matrix [B0;B1;...;Bp]
#   pf         cell array {x^0,x^1,...,x^p}
#   C          cell array {C1,C2,...,Cq}
#   CC         matrix [C1;C2;...;Cq]
#   f          cell array {f1,f2,...,fq}
#   iL         vector with indices of L-factors
#   L          cell array {L1,L2,...,Lq}
#   LL         matrix [L1,L2,...,Lq]
#   UU         matrix [U1,U2,...,Uq]
#   n          dimension of NLEP
#   p          order of polynomial part (length of B = p + 1)
#   q          length of f and C
#   r          sum of ranks of low rank matrices in C
#   Sigma      vector with target set (polygon)
#   leja       0: no leja, 1: leja in expansion phase, 2: always leja
#   nodes      [] or given nodes
#   Xi         vector with singularity set (discretized)
#   computeD   is true if generalized divided differences are explicitly used
#   tollin     tolerance for convergence of linearization
#   tolres     tolerance for residual
#   maxdgr     maximum degree of approximation
#   minit      minimum number of iterations after linearization is converged
#   maxit      maximum number of iterations
#   isfunm     is true if f are matrix functions
#   static     is true if static version is used
#   v0         starting vector
#   reuselu    positive integer for reuse of LU-factorizations of A(sigma)
#   funres     function handle for residual R(Lambda,X)
#   b          block size for pre-allocation
#   verbose    level of display [ {0} | 1 | 2 ]
function prepare_inputs(nep::NEP, Sigma::AbstractVector{Complex{T}}, Xi::AbstractVector{T}, options::Dict) where T<:Real
    iL = Vector{Int}(0)
    L = Vector{Matrix{T}}(0)
    LL = Matrix{T}(0, 0)
    UU = Matrix{T}(0, 0)
    #L = Vector{eltype(nep.spmf.A)}(0)
    #LL = similar(nep.spmf.A[1], 0, 0)
    #UU = similar(nep.spmf.A[1], 0, 0)
    p = -1
    q = 0
    r = 0
    spmf = isa(nep, AbstractSPMF)

    if !spmf
        BBCC = Matrix{T}(0, 0)
        #BBCC = isempty(nep.B) ? similar(nep.C[1].A, 0, 0) : similar(nep.B[1], 0, 0)
    else
        p = nep.p
        q = nep.q
        if q > 0
            # L and U factors of the low rank nonlinear part C
            L = nep.L
            if !isempty(nep.L[1])
                UU = hcat(nep.U...)::eltype(nep.U)
                r = nep.r
                iL = zeros(Int, r)
                c = 0
                for ii = 1:q
                    ri = size(nep.L[ii], 2)
                    iL[c+1:c+ri] = ii
                    c += ri
                end
                LL = hcat(nep.L...)::eltype(nep.L)
            end
        end

        BBCC = vcat(nep.spmf.A...)::eltype(nep.spmf.A)
    end

    # process the input Xi
    if isempty(Xi)
        Xi = [Inf]
    end

    # extract input options, with default values if missing
    verbose = get(options, "disp", 0)::Int
    maxdgr = get(options, "maxdgr", 100)::Int
    minit = get(options, "minit", 20)::Int
    maxit = get(options, "maxit", 200)::Int
    tolres = get(options, "tolres", 1e-10)::T
    tollin = get(options, "tollin", max(tolres/10, 100*eps()))::T
    v0 = get(options, "v0", randn(nep.n))::Vector{T}
    residual = [] # TODO
    funres = get(options, "funres", residual) # TODO: populate with default residual
    isfunm = get(options, "isfunm", true)::Bool
    static = get(options, "static", false)::Bool
    leja = get(options, "leja", 1)::Int
    nodes = get(options, "nodes", Vector{Complex{T}}(0))::Vector{Complex{T}}
    reuselu = get(options, "reuselu", 1)::Int
    b = get(options, "blksize", 20)::Int

    # extra defaults
    computeD = (nep.n <= 400) # for small problems, explicitly use generalized divided differences
    resfreq = 5

    if nep.n == 1
        maxdgr = maxit + 1;
    end

    return spmf,iL,L,LL,UU,p,q,r,Sigma,leja,nodes,
            Xi,tollin,tolres,maxdgr,minit,maxit,isfunm,static,v0,reuselu,
            funres,b,computeD,resfreq,verbose,BBCC
end

# scgendivdiffs: compute scalar generalized divided differences
#   sigma   discretization of target set
#   xi      discretization of singularity set
#   beta    scaling factors
function scgendivdiffs(sigma, xi, beta, maxdgr, isfunm, pff)
    sgdd = complex(zeros(length(pff), maxdgr+2))
    for ii = 1:length(pff)
        if isfunm
            sgdd[ii,:] = ratnewtoncoeffsm(pff[ii], sigma, xi, beta)
        else
            sgdd[ii,:] = map(m -> m[1], ratnewtoncoeffs(pff[ii], sigma, xi, beta))
        end
    end
    return sgdd
end

# constructD: Construct generalized divided difference
#   nb  number
function constructD(nb, L, n, p, q, r, BC, sgdd)
    if r == 0 || nb <= p
        D = spzeros(n, n)
        for ii = 1:(p+1+q)
            D += sgdd[ii,nb+1] * BC[ii]
        end
    else
        D = []
        for ii = 1:q
            d = sgdd[p+1+ii,nb+1] * L[ii]
            D = ii == 1 ? d : hcat(D, d)
        end
    end
    return D
end

# backslash: Backslash or left matrix divide
#   wc       continuation vector
function backslash(wc, nep, spmf, iL, LL, UU, p, q, r, reuselu, computeD, sigma, k, D, beta, N, xi, expand, kconv, BBCC, sgdd)
    n = nep.n
    shift = sigma[k+1]

    # construction of B*wc
    Bw = zeros(eltype(wc), size(wc))
    # first block (if low rank)
    if r > 0
        i0b = (p-1)*n + 1
        i0e = p*n
        if !spmf || computeD
            Bw[1:n] = -D[p+1] * wc[i0b:i0e] / beta[p+1]
        else
            Bw[1:n] = -sum(reshape(BBCC * wc[i0b:i0e], n, :) .* sgdd[:,p+1].', 2) / beta[p+1];
        end
    end
    # other blocks
    i0b = 1
    i0e = n
    for ii = 1:N
        # range of block i+1
        i1b = i0e + 1
        if r == 0 || ii < p
            i1e = i0e + n
        else
            i1e = i0e + r
        end
        # compute block i+1
        if r == 0 || ii != p
            Bw[i1b:i1e] = wc[i0b:i0e] + beta[ii+1]/xi[ii]*wc[i1b:i1e]
        else
            Bw[i1b:i1e] = UU' * wc[i0b:i0e] + beta[ii+1]/xi[ii]*wc[i1b:i1e]
        end
        # range of next block i
        i0b = i1b
        i0e = i1e
    end

    # construction of z0
    z = copy(Bw)
    i1b = n + 1
    if r == 0 || p > 1
        i1e = 2*n
    else
        i1e = n + r
    end
    nu = beta[2]*(1 - shift/xi[1])
    z[i1b:i1e] = 1/nu * z[i1b:i1e]
    for ii = 1:N
        # range of block i+2
        i2b = i1e + 1
        if r == 0 || ii < p-1
            i2e = i1e + n
        else
            i2e = i1e + r
        end
        # add extra term to z0
        if !spmf || computeD
            if r == 0 || ii != p
                z[1:n] -= D[ii+1] * z[i1b:i1e]
            end
        else
            if r == 0 || ii < p
                z[1:n] -= sum(reshape(BBCC * z[i1b:i1e], n, :) .* sgdd[:,ii+1].', 2)
            elseif ii > p
                dd = sgdd[p+2:end,ii+1]
                z[1:n] -= LL*(z[i1b:i1e] .* dd[iL])
            end
        end
        # update block i+2
        if ii < N
            mu = shift - sigma[ii+1]
            nu = beta[ii+2] * (1 - shift/xi[ii+1])
            if r == 0 || ii != p-1
                z[i2b:i2e] = 1/nu * z[i2b:i2e] + mu/nu * z[i1b:i1e]
            else # i == p-1
                z[i2b:i2e] = 1/nu * z[i2b:i2e] + mu/nu * UU'*z[i1b:i1e]
            end
        end
        # range of next block i+1
        i1b = i2b
        i1e = i2e
    end

    # solving Alam x0 = z0
    w = zeros(eltype(wc), size(wc))
    if ((!expand || k > kconv) && reuselu == 1) || reuselu == 2
        w[1:n] = lusolve(λ -> compute_Mder(nep, λ), shift, z[1:n]/beta[1])
    else
        w[1:n] = compute_Mder(nep, shift) \ (z[1:n]/beta[1])
    end

    # substitutions x[i+1] = mu/nu*x[i] + 1/nu*Bw[i+1]
    i0b = 1
    i0e = n
    for ii = 1:N
        # range of block i+1
        i1b = i0e + 1
        if r == 0 || ii < p
            i1e = i0e + n
        else
            i1e = i0e + r
        end
        # compute block i+1
        mu = shift - sigma[ii]
        nu = beta[ii+1] * (1 - shift/xi[ii])
        if r == 0 || ii != p
            w[i1b:i1e] = mu/nu * w[i0b:i0e] + 1/nu * Bw[i1b:i1e]
        else
            w[i1b:i1e] = mu/nu * UU'*w[i0b:i0e] + 1/nu * Bw[i1b:i1e]
        end
        # range of next block i
        i0b = i1b
        i0e = i1e
    end

    return w
end

# in_sigma: True for points inside Sigma
#   z      (complex) points
function in_sigma(z::AbstractVector{Complex{T}}, Sigma::AbstractVector{Complex{T}}, tolres::T) where T<:Real
    if length(Sigma) == 2 && isreal(Sigma)
        realSigma = real([Sigma[1]; Sigma[1]; Sigma[2]; Sigma[2]])
        imagSigma = [-tolres; tolres; tolres; -tolres]
    else
        realSigma = real(Sigma)
        imagSigma = imag(Sigma)
    end
    return map(p -> inpolygon(real(p), imag(p), realSigma, imagSigma), z)
end

#=
# residual: Residual of the NLEP for given eigenvalues and eigenvectors
#   Lambda  eigenvalues
#   X       eigenvectors (normalized)
    function residual(Lambda, X, funA, Ahandle, BB, pf, CC, f, p, q)
        R = zeros(length(Lambda), 1)
        if Ahandle
            for ii = 1:length(Lambda)
                R[ii] = norm(funA(Lambda[ii]) * X[:,ii])
            end
        else
            if p < 0
                BCX = reshape(CC*X, [], q, size(X,2))
            else
                BCX = reshape([BB*X; CC*X], [], p+q+1, size(X,2))
            end
            FL = cell2mat(cellfun(@(x) arrayfun(x,Lambda),[pf[:]; f[:]]', 'UniformOutput',0))
            RR = sum(bsxfun(@times, BCX, reshape(FL.', 1, size(BCX, 2), [])), 2)
            for ii = 1:length(Lambda)
                R[ii] = norm(RR[:, :, ii])
            end
        end
        return R
    end
=#

function resize_matrix(A, rows, cols)
    resized = zeros(eltype(A), rows, cols)
    resized[1:size(A, 1), 1:size(A, 2)] = A
    return resized
end
