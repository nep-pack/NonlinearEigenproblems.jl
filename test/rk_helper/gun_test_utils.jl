using NonlinearEigenproblems
using NonlinearEigenproblems.RKHelper
using Random
using LinearAlgebra

function gun_init()
    nep = gun_nep()

    gam = 300^2 - 200^2
    mu = 250^2
    sigma2 = 108.8774
    xmin = gam*(-1) + mu
    xmax = gam*1 + mu

    # define target set Σ
    npts = 1000
    halfcircle = xmin .+ (xmax-xmin) * (cis.(range(0, stop = pi, length = round(Int, pi/2*npts) + 2)) / 2 .+ .5)
    Σ = [halfcircle; xmin]

    # sequence of interpolation nodes
    Z = [2/3, (1+im)/3, 0, (-1+im)/3, -2/3]
    nodes = gam*Z .+ mu

    # define the set of pole candidates
    Ξ = -10 .^ range(-8, stop = 8, length = 10000) .+ sigma2^2

    # options
    Random.seed!(1)
    v = randn(size(nep, 1)) .+ 0im

    funres = (λ, v) -> gun_residual(λ, v, nep.nep1.A..., nep.nep2.spmf.A...)

    return nep, Σ, Ξ, v, nodes, funres
end

function gun_nep()
    nep = nep_gallery("nlevp_native_gun")
    K, M = get_Av(nep.nep1);
    W1, W2 = get_Av(nep.nep2);
    c1 = LowRankMatrixAndFunction(W1, get_fv(nep.nep2)[1])
    c2 = LowRankMatrixAndFunction(W2, get_fv(nep.nep2)[2])
    return SumNEP(PEP([K, M]), LowRankFactorizedNEP([c1, c2]))
end

function gun_residual(λ, v, K, M, W1, W2)
    # constants
    sigma1 = 0
    sigma2 = 108.8774

    nK = 1.474544889815002e+05   # opnorm(K, 1)
    nM = 2.726114618171165e-02   # opnorm(M, 1)
    nW1 = 2.328612251920476e+00  # opnorm(W1, 1)
    nW2 = 3.793375498194695e+00  # opnorm(W2, 1)

    # Denominator
    den = nK + abs(λ) * nM + sqrt(abs(λ-sigma1^2)) * nW1 + sqrt(abs(λ-sigma2^2)) * nW2

    # 2-norm of A(lambda)*x
    norm((K + M*λ + W1*im*sqrt(λ) + W2*im*sqrt(λ - sigma2^2)) * v) / den
end
