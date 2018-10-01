using NonlinearEigenproblems.RKHelper
using LinearAlgebra
using Random

export cork, CorkSolutionDetails

"""
    cork(nep::NEP)

Find a few eigenvalues and eigenvectors of a linearization pencil.

# Arguments
- `nep`: An instance of a nonlinear eigenvalue problem.
- `L`: Linearization pencil. See Linearization.
- `σ`: Scalar for finding the eigenvalues closest to.
- `Σ`: A vector containing the points of a polygonal target set in the complex plane.
- `m`: Maximum number of ritz values.
- `p`: Number of restarted ritz values.
- `shifts`: Cyclically repeated shifts for the rational Krylov process.
- `maxrest`: Maximum number or restarts.
- `displaylevel`: Level of display (0, 1, 2).
- `tol`: Tolerance for residual.
- `tolrnk`: Tolerance for truncating rank.
- `tollck`: Tolerance for locking.
- `v`: Starting vector (sized any multiple of NEP size).
- `errmeasure`: Function for error measure (residual norm). Called with arguments (λ,v).
- `reusefact`: Whether to reuse matrix factorizations.
- `return_details`: Whether to return solution details (see NleigsSolutionDetails).
- `k`: Desired number of eigenvalues.
   cork(L) returns a vector lambda of the 6 smallest magnitude
   cork(L,k) returns the k largest magnitude eigenvalues.
   cork(L,k,sigma) returns the k eigenvalues closest to the real or complex scalar sigma.
   cork(L,k,target) returns the eigenvalues inside the target set target.Sigma
   and closest to the shift target.sigma:
     target.Sigma:  vector containing the points of a polygonal target set in
                    the complex plane
     target.shift:  scalar shift
   k is a guess of the number of eigenvalues inside the target set.

# Return values
- `λ`: Vector of eigenvalues of the nonlinear eigenvalue problem NLEP inside the target set Σ.
- `X`: Corresponding matrix of eigenvectors.
- `res`: Corresponding residuals.
- `flag`: Convergence flag. True if all the eigenvalues converged.
- `details`: Solution details, if requested (see CorkSolutionDetails).

# References
- R. Van Beeumen, K. Meerbergen, and W. Michiels. Compact rational Krylov
  methods for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 36(2),
  pp. 820-838, 2015.
- [CORK Matlab toolbox](http://twr.cs.kuleuven.be/research/software/nleps/cork.php)
"""
function cork(
    ::Type{T},
    nep::NEP,
    L::Linearization,
    σ::CT = 0,
    Σ::AbstractVector{CT} = Vector{CT}();
    k::Int = 6,
    m::Int = max(3*k, 20),
    p::Int = max(2*k, 10),
    shifts::AbstractVector{CT} = [σ],
    maxrest::Int = 50,
    displaylevel::Int = 0,
    linsolvercreator::Function = default_linsolvercreator,
    tolres::T = 100*eps(T),
    tolrnk::Vector{T} = Vector{T}(),
    tollck::T = eps(T),
    v0::Vector{CT} = CT.(randn(T, size(nep, 1))),
    errmeasure::Function = default_errmeasure(nep::NEP),
    return_details::Bool = false) where {T<:Real, CT<:Complex{T}}

    n = size(nep, 1)
    d = length(L.A)
    reuselu = false
    maxreorth = 2

    nb = ceil(Int, (m + maxrest*(m-p) + 1) / length(shifts))
    shifts = repeat(reshape(shifts, :, 1), nb, 1)

    # initialization
    Q,U,r,j,H,K,X,lambda,res,Lam,Res,J,R,LU,SHIFTS,NNZ = initialize(CT, v0, L, n, m, d, p, maxrest, return_details)

    # CORK loop
    i = 1
    nbrest = 0
    l = 0
    count = 0
    flag = true

    while i <= m + maxrest*(m-p)
        r = corkstep(shifts[i], r, j, i, L, Q, U, d, NNZ, displaylevel)

        return ([],[],[],true,CorkSolutionDetails{T,CT}())

        # compute Ritz pairs and residuals
        flag = ritzpairs(r, j, i, l)

        # check for convergence
        !flag && break

        # restart
        if j == m
            if nbrest < maxrest
                nbrest += 1
                r,j,l = implicitrestart(r, nbrest)
            else
                break
            end
        end

        i += 1
        j += 1
    end

    outputs(return_details, flag, i)
end

function initialize(CT, v0, L, n, m, d, p, maxrest, return_details)
    # matrix Q -- start vector can be sized any multiple of 'n'
    Q,U0 = qr(reshape(normalize(v0), n, :))
    Q = Matrix(Q) # Thin QR-factorization
    r = size(Q, 2)

    # tensor U
    U = zeros(CT, r, m+1, d)
    U[1:r,1,1:size(U0,1)] = U0
    j = 1

    # Hessenberg matrices H and K
    H = zeros(m+1, m)
    K = zeros(m+1, m)

    # eigenpairs and residuals
    X = fill(NaN, n, m)
    lambda = fill(NaN, m)
    res = fill(NaN, m)

    # info variables
    if return_details
        Lam = fill(NaN, m, m + maxrest*(m-p))
        Res = fill(NaN, m, m + maxrest*(m-p))
    else
        Lam = []
        Res = []
    end
    J = [1; zeros(m + maxrest*(m-p))]
    R = [r; zeros(m + maxrest*(m-p))]

    # lu variables
    LU = []
    SHIFTS = []
    anynnz(X) = nnz(X) > 0
    NNZ = (A = map(anynnz, L.A), B = map(anynnz, L.B))

    return (Q,U,r,j,H,K,X,lambda,res,Lam,Res,J,R,LU,SHIFTS,NNZ)
end

function corkstep(shift, r, j, i, L, Q, U, d, NNZ, displaylevel)
    displaylevel > 0 && println("iteration $i")
    eta = 1/sqrt(2) # reorthogonalization constant

    # compute next vector
    v = Q[:, 1:r] * reshape(U[1:r, j, 1:d], r, :)
    v1,MsN1inv0,MsN1invN = backslash(shift, v, L, d, NNZ)

    # level 1 orthogonalization
    q = v1
    delta = Inf
    nborth = 0
    while norm(q) < eta*delta && nborth <= maxreorth
        delta = norm(q)
        if nborth == 0
            u1 = Q[:,1:r]'*v1
            q = q - Q[:,1:r]*u1
        else
            q = q - Q[:,1:r] * (Q[:,1:r]'*q)
        end
        nborth += 1
    end
    delta = norm(q)
    # update Q
    if delta > eps(T)
        rnew = r + 1
        q = q/delta
        Q[:,rnew] = q
    else
        rnew = r
    end

    # level 2 orthogonalization
    if rnew > r
        U[rnew,:,:] = 0
        u1 = [u1; q'*v1]
    end
    u2 = [] #TODO: reshape(U(1:rnew,j,1:d),rnew,[])*transpose(MsN1invN) - bsxfun(@times,transpose(MsN1inv0),u1)
    u = [u1; reshape(u2, :, 1)]
    delta = Inf
    nborth = 0
    while norm(u) < eta*delta && nborth <= maxreorth
        delta = norm(u)
        h = reshape(conj(permute(U[1:rnew,1:j,1:d], [2,1,3])), :, d*rnew) * u
        u = u - reshape(permute(U[1:rnew,1:j,1:d], [1,3,2]), :, j) * h
        H[1:j,j] = H[1:j,j] + h
        nborth += 1
    end
    K[1:j,j] = shift*H[1:j,j] + flipud(eye(j,1))
    H[j+1,j] = norm(u)
    K[j+1,j] = shift*H[j+1,j]
    # update U
    u = u/H[j+1,j]
    U[1:rnew,j+1,1:d] = reshape(u, :, d)

    # dimensions
    J[i+1] = j+1
    R[i+1] = rnew

    return rnew
end

function backslash(shift, y, L, d, NNZ)
    # check for new shift
    newlu = true
    # TODO: the below will be replaced by LinSolverCache
#=    if !reuselu
        newlu = true
    else
        newlu = false
        if isempty(SHIFTS)
            idx = 1;
            SHIFTS[idx] = shift
            newlu = true
        else
            idx = find(SHIFTS == shift, 1, "first")
            if isempty(idx)
                idx = length(SHIFTS) + 1
                SHIFTS[idx] = shift
                newlu = true
            end
        end
    end=#

    # initialize
    m0 = L.M[1:d-1,1]
    M1 = L.M[1:d-1,2:d]
    n0 = L.N[1:d-1,1]
    N1 = L.N[1:d-1,2:d]

    # compute kronecker coefficients
    MsN1 = M1 - shift * N1
    msn0 = m0 - shift * n0
    MsN1inv0 = MsN1 \ msn0
    MsN1invN = MsN1 \ L.N[1:d-1,1:d]

    # build At
    if newlu
        if !ismissing(L.Pfun)
            At = L.Pfun(shift)
        else
            At = L.A[1] - shift*L.B[1]
            for ii = 2:d
                At = At - MsN1inv0[ii-1]*(L.A[ii] - shift*L.B[ii])
            end
        end
    end

    # intermediate x used in yt
    x = y * transpose(MsN1invN)

    # build yt
    yt = L.B[1] * y[:,1]
    for ii = 2:d
        if NNZ.A[ii] && NNZ.B[ii]
            # experimental:
#            corr = L.A{ii}*x(:,ii-1) - L.B{ii}*(shift*x(:,ii-1)+y(:,ii));
#            norm(corr) < tolres && break
#            yt -= corr
            yt -= L.A[ii] * x[:,ii-1] - L.B[ii] * (shift*x[:,ii-1] + y[:,ii])
        elseif NNZ.A[ii]
            yt -= L.A[ii] * x[:,ii-1]
        else
            yt = yt + L.B[ii] * (shift*x[:,ii-1] + y[:,ii])
        end
    end

    # compute x1 = At\yt
    # TODO: use linsolvercreator
#    if !reuselu
        x1 = At \ yt
        #=
    else
        if newlu
            if issparse(At)
                AL,AU,AP,AQ,AR = lu(At)
                x1 = AQ * (AU \ (AL \ (AP * (AR \ yt))));
                % improve accuracy
                resid = yt - At*x1;
                err = AQ * (AU \ (AL \ (AP * (AR \ resid))));
                x1 = x1 + err;
                LU{idx} = struct('A',At,'L',AL,'U',AU,'P',AP,'Q',AQ,'R',AR);
            else
                [AL,AU,Ap] = lu(At,'vector');
                x1 = AU \ (AL \ yt(Ap,:));
                % improve accuracy
                resid = yt - At*x1;
                err = AU \ (AL \ resid(Ap,:));
                x1 = x1 + err;
                LU{idx} = struct('A',At,'L',AL,'U',AU,'p',Ap);
            end
        else
            if isfield(LU{idx},'R')
                x1 = LU{idx}.Q * (LU{idx}.U \ (LU{idx}.L \ ...
                    (LU{idx}.P * (LU{idx}.R \ yt))));
                % impove accuracy
                resid = yt - LU{idx}.A*x1;
                err = LU{idx}.Q * (LU{idx}.U \ (LU{idx}.L \ ...
                    (LU{idx}.P * (LU{idx}.R \ resid))));
                x1 = x1 + err;
            else
                x1 = LU{idx}.U \ (LU{idx}.L \ yt(LU{idx}.p,:));
                % improve accuracy
                resid = yt - LU{idx}.A*x1;
                err = LU{idx}.U \ (LU{idx}.L \ resid(LU{idx}.p,:));
                x1 = x1 + err;
            end
        end
    end
    =#

    return (x1,MsN1inv0,MsN1invN)
end


function ritzpairs(r,j,i,l)
    # QZ factorization
    qzreal = isreal(σ) && isreal(K) && isreal(H)
    if qzreal
        G,F,_,_,S,_ = qz(K[1:j,1:j], H[1:j,1:j], "real")
    else
        G,F,_,_,S,_ = qz(K[1:j,1:j], H[1:j,1:j])
    end

    # ritz pairs
    ordlam = ordeig(G, F)
    lambda[l+1:j] = ordlam(l+1:j)
    idx,in = sortlam(lambda[1:j])
    if qzreal
        ii = l+1
        while ii <= j
            if isreal(lambda(ii))
                S[:,ii] = S[:,ii] / norm(S[:,ii])
                ii += 1
            else
                Sr = S[:,ii]
                Si = S[:,ii+1]
                S[:,ii] = (Sr + 1i*Si) / norm(Sr + 1i*Si)
                S[:,ii+1] = (Sr - 1i*Si) / norm(Sr - 1i*Si)
                ii += 2
            end
        end
    else
        for ii = l+1:j
            S[:,ii] = S[:,ii] / norm(S[:,ii])
        end
    end
    X[:,l+1:j] = Q[:,1:r] * (U[1:r,1:j+1,1] * (H[1:j+1,1:j] * S[:,l+1:j]))
    for ii = l+1:j
        X[:,ii] = X[:,ii] / norm(X[:,ii])
    end

    # residuals
    if isempty(funres)
        res[l+1:j] = abs((K[j+1,j] - lambda[l+1:j]*H[j+1,j]) .* transpose(S[j,l+1:j]))
    else
        res[l+1:end] = NaN
        in2 = in
        in2[1:l] = false
        res[in2] = funres(lambda[in2], X[:,in2])
    end

    # info
    if info
        Lam[1:j,i] = lambda[1:j]
        Res[1:j,i] = res[1:j]
    end

    # check for convergence
    if isempty(Σ)
        if length(lambda[1:j]) < k
            count = 0
        else
            if all(res(idx[1:k]) <= tolres)
                count += 1
            else
                count = 0
            end
        end
    else
        if nbrest == 0
            count = 0
        else
            if all(res[in] <= tolres)
                count += 1
            else
                count = 0
            end
        end
    end
    flag = count < 5

    if displaylevel > 0
        if isempty(Σ)
            @printf(" %i < %.2e\n", sum(res[idx[1:min(j,k)]] <= tolres), tolres)
        else
            @printf(" %i (of %i) < %.2e\n", sum(res[in] <= tolres), sum(in), tolres)
        end
        if !flag
            if isempty(Σ)
                println(" ==> $k converged Ritz pairs")
            else
                println(" ==> all Ritz pairs inside target set converged")
            end
        end
    end

    # convergence
    if !flag
        if isempty(Σ)
            lambda = lambda[idx[1:k]]
            X = X[:,idx[1:k]]
            res = res[idx[1:k]]
        else
            lambda = lambda[in]
            X = X[:,in]
            res = res[in]
        end
    end

    return flag
end

function sortlam(lambda)
    if isempty(Σ)
        _,idx = sort(abs(lambda - σ))
        in = fill(true, length(lambda))
    else
        in = inpolygon(real(lambda), imag(lambda), real(Σ), imag(Σ))
        Ilambda = 1:length(lambda)
        Ilamin = Ilambda[in]
        Ilamout = Ilambda[!in]
        _,Iin = sort(abs(lambda[in] - σ))
        _,Iout = sort(abs(lambda[!in] - σ))
        idx = [Ilamin[Iin]; Ilamout[Iout]]
    end
    return (idx,in)
end

function implicitrestart(r, nbrest)
    # QZ factorization
    qzreal = isreal(σ) && isreal(K) && isreal(H)
    if qzreal
        G,F,Y,Z = qz(K[1:m,:], H[1:m,:], "real")
    else
        G,F,Y,Z = qz(K[1:m,:], H[1:m,:])
    end

    # select ritz values and reorder QZ factorization
    Ilambda = sortlam(lambda)
    if qzreal && !isreal(lambda[Ilambda[p]]) && imag(lambda[Ilambda[p]]) == -imag(lambda[Ilambda[p+1]])
        # increment p in case of complex conjugated ritz value
        j = p + 1
    else
        j = p
    end
    selres = zeros(m, 1)
    _,Ires = sort(res)
    selres[Ires] = m:-1:1
    selres[res == 0] = 2*m
    select = zeros(m, 1)
    if qzreal
        ii = 0
        while ii < j
            if isreal(lambda[Ilambda[ii+1]])
                select[Ilambda[ii+1]] = selres[Ilambda[ii+1]]
                ii += 1
            else
                select[Ilambda[ii+1:ii+2]] = selres[Ilambda[ii+1]]
                ii += 2
            end
        end
    else
        select[Ilambda[1:j]] = selres[Ilambda[1:j]]
    end
    G,F,Y,Z = ordqz(G, F, Y, Z, select)

    # new Hessenberg matrices H and K
    Hj1 = H[m+1,:] * Z[:,1:j]
    Kj1 = K[m+1,:] * Z[:,1:j]
    # locking
    lnew = find(abs(Hj1) > tollck, 1, "first") - 1
    if lnew > 0
        _,Ilock = sort(select, "descend")
        if qzreal && !isreal(lambda[Ilock[lnew]]) && imag(lambda[Ilock[lnew]]) > 0
            # decrement lnew in case of missing complex conjugated value
            lnew -= 1
        end
        lambda[1:lnew] = lambda[Ilock[1:lnew]]
        X[:,1:lnew] = X[:,Ilock[1:lnew]]
        res[1:lnew] = res[Ilock[1:lnew]]
        Hj1[1:lnew] = 0
        Kj1[1:lnew] = 0
    end
    H = zeros(m+1,m)
    K = zeros(m+1,m)
    H[1:j+1,1:j] = [F[1:j,1:j]; Hj1]
    K[1:j+1,1:j] = [G[1:j,1:j]; Kj1]

    # new tensor U
    U = permute(reshape([reshape(permute(U[:,1:m,:],[1,3,2]),[],m) * Y(1:j,:)',reshape(U[:,m+1,:],[],1)],r,[],j+1),[1,3,2]);
    UU,US,_ = svd(reshape(U,r,[]), "econ")
    Us = diag(US)
    if isempty(tolrnk)
        rnew = sum(Us > max(r,(j+1)*d) * eps(max(Us))) # see rank.m
    else
        rnew = sum(Us > tolrnk)
    end
    U = reshape(UU[:,1:rnew]'*reshape(U[1:r,1:j+1,1:d],r,[]),rnew,[],d)

    # new matrix Q
    Q = Q*UU[:,1:rnew]

    displaylevel > 0 && println(" ==> restart $nbrest: r = $rnew, j = $j, l = $lnew")

    return (rnew,j,lnew)
end

function outputs(return_details, flag, i)
    # no convergence
    if flag
        idx,in = sortlam(lambda)
        if isempty(Σ)
            lambda = lambda[idx[1:k]]
            X = X[:,idx[1:k]]
            res = res[idx[1:k]]
        else
            lambda = lambda[in]
            X = X[:,in]
            res = res[in]
        end
        # number of converged eigenvalues
        nb = sum(res <= tolres)
        if isempty(Σ)
            lambda = [lambda[res <= tolres]; fill(NaN, k-nb)]
            X = [X[:,res <= tolres], fill(NaN, n, k-nb)]
            res = [res[res <= tolres]; fill(NaN, k-nb)]
        end
        if nb == 0
            @warn "CORK: None of the $k requested eigenvalues converged."
        else
            @warn "CORK: Only $nb of the $k requested eigenvalues converged."
        end
    end

    details = return_details ?
        CorkSolutionDetails(Lam[:,1:i], Res[:,1:i], J[1:i+1], R[1:i+1], m, p) :
        CorkSolutionDetails{T,CT}()

    return lambda, X, res, flag, details
end

struct CorkSolutionDetails{T<:Real, CT<:Complex{T}}
    "matrix of Ritz values in each iteration"
    Lam::Matrix{CT}

    "matrix of residuals in each iteraion"
    Res::Matrix{T}

    "vector with dimension of Krylov subspace in each iteration"
    J::Vector{Int}

    "vector with rank of subspace Q in each iteration"
    R::Vector{Int}

    "maximum number of ritz values"
    m::Int

    "number of restarted ritz values"
    p::Int
end

CorkSolutionDetails{T,CT}() where {T<:Real, CT<:Complex{T}} = CorkSolutionDetails(
    Matrix{CT}(undef, 0, 0), Matrix{T}(undef, 0, 0), Vector{Int}(), Vector{Int}(), 0, 0)
