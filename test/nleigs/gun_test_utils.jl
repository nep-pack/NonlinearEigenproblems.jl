function gun_init()
    nep = gun_nlep()

    gam = 300^2 - 200^2
    mu = 250^2
    sigma2 = 108.8774
    xmin = gam*(-1) + mu
    xmax = gam*1 + mu

    # define target set Sigma
    npts = 1000
    halfcircle = xmin + (xmax-xmin) * (cis.(linspace(0, pi, round(pi/2*npts) + 2)) / 2 + .5)
    Sigma = [halfcircle; xmin]

    # sequence of interpolation nodes
    Z = [2/3, (1+im)/3, 0, (-1+im)/3, -2/3].'
    nodes = gam*Z + mu

    # define the set of pole candidates
    Xi = -logspace(-8, 8, 10000) + sigma2^2

    # options
    srand(1)
    v0 = randn(nep.n)

    funres = (Lam, X) -> gun_residual(Lam, X, nep.B[1], nep.B[2], nep.C[1].A, nep.C[2].A)

    return nep, Sigma, Xi, v0, nodes, funres
end

function gun_nlep()
    nep = nep_gallery("nlevp_native_gun")

    K = nep.A[1]
    M = nep.A[2]
    W1 = nep.A[3]
    W2 = nep.A[4]

    sigma2 = 108.8774

    C = [W1, W2] # nonlinear matrices

    # compute svd decomposition of C[1]
    W1 = C[1][9722:end, 9722:end]
    #lu = lufact(W1) # TODO use this
    L, U = lu(full(W1))
    L1, U1 = compactlu(sparse(L), sparse(U))
    L1a = spzeros(9956, size(L1, 2))
    L1a[9722:end, :] = L1
    U1a = spzeros(size(U1, 1), 9956)
    U1a[:, 9722:end] = U1
    U1a = U1a'

    # compute svd decomposition of C[2]
    W2 = C[2][1:479, 1:479]
    #lu = lufact(W2) # TODO use this
    L, U = lu(full(W2))
    L2, U2 = compactlu(sparse(L), sparse(U))
    L2a = spzeros(9956, size(L2, 2));
    L2a[1:479, :] = L2
    U2a = spzeros(size(U2, 1), 9956)
    U2a[:, 1:479] = U2
    U2a = U2a'

    # nonlinear functions
    f = [nep.fi[3], nep.fi[4]]

    # finally assemble nep instance
    c1 = SPMFLowRankMatrix(C[1], L1a, U1a, f[1])
    c2 = SPMFLowRankMatrix(C[2], L2a, U2a, f[2])
    return SPMFLowRankNEP(9956, [K, -M], [c1, c2])
end

function compactlu(L, U)
    n = size(L, 1)

    select = map(i -> nnz(L[i:n, i]) > 1 || nnz(U[i, i:n]) > 0, 1:n)

    Lc = L[:,select]
    Uc = U[select,:]

    return Lc, Uc
end

function gun_residual(Lambda, X, K, M, W1, W2)
    # constants
    sigma1 = 0
    sigma2 = 108.8774

    nK = 1.474544889815002e+05   # norm(K, 1)
    nM = 2.726114618171165e-02   # norm(M, 1)
    nW1 = 2.328612251920476e+00  # norm(W1, 1)
    nW2 = 3.793375498194695e+00  # norm(W2, 1)

    # Denominator
    Den = nK + abs.(Lambda) * nM + sqrt.(abs.(Lambda-sigma1^2)) * nW1 + sqrt.(abs.(Lambda-sigma2^2)) * nW2

    # 2-norm of A(lambda)*x
    R = map(i -> norm((K + M*Lambda[i] + W1*im*sqrt(Lambda[i]) + W2*im*sqrt(Lambda[i] - sigma2^2)) * X[:,i]) / Den[i], 1:length(Lambda))

    return R
end
