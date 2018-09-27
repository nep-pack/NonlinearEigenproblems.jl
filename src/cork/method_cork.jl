using NonlinearEigenproblems

include("../nleigs/method_nleigs.jl")

function cork()
    nep = nep_gallery("nlevp_native_gun")
    K, M = get_Av(nep.nep1);
    W1, W2 = get_Av(nep.nep2);
    c1 = LowRankMatrixAndFunction(W1, get_fv(nep.nep2)[1])
    c2 = LowRankMatrixAndFunction(W2, get_fv(nep.nep2)[2])
    nep = SumNEP(PEP([K, M]), LowRankFactorizedNEP([c1, c2]))

    T = Float64
    P = get_nleigs_nep(T, nep)

    # Gun-specific setup
    gam = 300^2 - 200^2
    mu = 250^2
    sigma2 = 108.8774
    xmin = gam*(-1) + mu
    xmax = gam*1 + mu

    # define target set Σ
    npts = 1000
    halfcircle = xmin .+ (xmax-xmin) * (cis.(range(0, stop = pi, length = round(Int, pi/2*npts) + 2)) / 2 .+ .5)
    Σ = [halfcircle; xmin]

    # define the set of pole candidates
    Ξ = -10 .^ range(-8, stop = 8, length = 10000) .+ sigma2^2

    L = linearize(T, P, Σ, Ξ, 100, 1e-10)
end

cork()
