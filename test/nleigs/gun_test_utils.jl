function gun_init()
    nep = nep_gallery("nlevp_gun_pnep")

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
    Z = [2/3, (1+im)/3, 0, (-1+im)/3, -2/3]
    nodes = gam*Z + mu

    # define the set of pole candidates
    Xi = -logspace(-8, 8, 10000) + sigma2^2

    # options
    srand(1)
    v0 = randn(nep.n)

    funres = (Lam, X) -> gun_residual(Lam, X, nep.spmf.A...)

    return nep, Sigma, Xi, v0, nodes, funres
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
