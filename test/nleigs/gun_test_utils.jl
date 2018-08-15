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
    v = randn(nep.n)

    funres = (λ, v) -> gun_residual(λ, v, nep.spmf.A...)

    return nep, Sigma, Xi, v, nodes, funres
end

function gun_residual(λ, v, K, M, W1, W2)
    # constants
    sigma1 = 0
    sigma2 = 108.8774

    nK = 1.474544889815002e+05   # norm(K, 1)
    nM = 2.726114618171165e-02   # norm(M, 1)
    nW1 = 2.328612251920476e+00  # norm(W1, 1)
    nW2 = 3.793375498194695e+00  # norm(W2, 1)

    # Denominator
    den = nK + abs(λ) * nM + sqrt(abs(λ-sigma1^2)) * nW1 + sqrt(abs(λ-sigma2^2)) * nW2

    # 2-norm of A(lambda)*x
    norm((K + M*λ + W1*im*sqrt(λ) + W2*im*sqrt(λ - sigma2^2)) * v) / den
end
