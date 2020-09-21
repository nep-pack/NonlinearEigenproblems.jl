using NonlinearEigenproblems
using NonlinearEigenproblems.RKHelper
using LinearAlgebra
using SparseArrays
using Random

function particle_init(interval)
    nep, brpts, U0 = particle_nep(interval)

    # interval
    sep = 1e-4
    if interval == 1
        xmin = -U0
        xmax = brpts[interval] - sep
        Ξ = 10 .^ range(-6, stop = 6, length = 10000) .+ brpts[interval]
    elseif interval > 1
        xmin = brpts[interval-1] + sep
        xmax = brpts[interval] - sep
        Ξ1 = -10 .^ range(-6, stop = 6, length = 5000) .+ brpts[interval-1]
        Ξ2 = 10 .^ range(-6, stop = 6, length = 5000) .+ brpts[interval]
        Ξ = [Ξ1; Ξ2]
    else
        error("Invalid interval: $interval")
    end

    # define target set Σ
    Σ = [xmin + 0im, xmax + 0im]

    # options
    temp_pep=nep_gallery("pep0",200); v=complex(vcat(vec(temp_pep.A[1][:,1:81]),temp_pep.A[1][1:81,82])) # Ugly but gives stability over versions. Stable "random" vector v.
    nodes = range(xmin + 0im, stop = xmax + 0im, length = 11)
    nodes = collect(nodes[2:2:end])

    return nep, Σ, Ξ, v, nodes, xmin, xmax
end

function particle_nep(interval)
    # Original MATLAB version by William Vandenberghe

    # initialization
    # constants
    meter = 1 / 5.2917725e-11
    nm = 1e-9 * meter
    eV = 1 / 13.6

    # parameters
    xmax = 5
    zmax = 2

    xstep = 0.05
    zstep = 0.05

    x_x = (-xmax:xstep:xmax) * nm
    z_z = (-zmax:zstep:zmax) * nm

    nx = length(x_x)
    nz = length(z_z)

    dx = minimum(diff(x_x))
    dz = minimum(diff(z_z))

    x = kron(x_x, z_z.^0)
    z = kron(x_x.^0, z_z)

    w1 = 1 * nm
    w2 = 1.1 * nm
    l = 4 * nm
    U0 = 3 * eV

    # potential U
    U = zero(x)
    U[abs.(z) .< w1] .= -U0
    U[(abs.(z) .< w2) .& (abs.(x) .< l/2)] .= -U0

    # m
    m = 0.2

    # branch points
    n = nx * nz
    Dxx_x = SymTridiagonal(ones(nx) * -2 / dx^2, ones(nx-1) / dx^2)
    Dzz_z = SymTridiagonal(ones(nz) * -2 / dx^2, ones(nz-1) / dz^2)
    H_L = -1/m*Dzz_z + diagm(0 => U[1:nz])
    H_R = -1/m*Dzz_z + diagm(0 => U[end-nz.+(1:nz)])

    if opnorm(H_L - H_R) == 0
        # symmetric potential
        D,V = eigen(H_L)
        i = sortperm(D)
        d = D[i]
        V = V[:,i]
        p = zeros(nz) # 0 (left+right)
    else
        # asymmetric potential
        D_L,V_L = eigen(H_L)
        D_R,V_R = eigen(H_R)
        V = [V_L; V_R]
        p = [-ones(nz); ones(nz)]
        DD = [D_L; D_R]
        i = sortperm(DD)
        d = DD[i]
        V = V[:,i]
        p = p[i] # -1 (left), 1 (right)
    end

    # matrix H
    H = -1/m * (kron(sparse(Dxx_x), sparse(1.0I, nz, nz)) + kron(sparse(1.0I, nx, nx), sparse(Dzz_z))) + sparse(Diagonal(U))

    # branch points
    brpts = unique(d)

    # polynomial matrices
    B = [H, sparse(-1.0I, n, n)]

    # nonlinear matrices
    SL = Vector(undef, length(brpts))
    if p[1] < 0
        SL[1] = [V[:,1]; spzeros(n-nz, 1)]
    elseif p[1] == 0
        SL[1] = hcat([V[:,1]; spzeros(n-nz, 1)], [spzeros(n-nz, 1); V[:,1]])
    else
        SL[1] = [spzeros(n-nz, 1); V[:,1]]
    end
    c = 1;
    for j = 2:length(d)
        if d[j-1] == d[j]
            if p[j] < 0
                SL[c] = hcat(SL[c], [V[:,j]; spzeros(n-nz, 1)])
            elseif p[j] == 0
                SL[c] = hcat(SL[c], [V[:,j]; spzeros(n-nz, 1)], [spzeros(n-nz, 1); V[:,j]])
            else
                SL[c] = hcat(SL[c], [spzeros(n-nz, 1); V[:,j]])
            end
        else
            c += 1
            if p[j] < 0
                SL[c] = [V[:,j]; spzeros(n-nz, 1)]
            elseif p[j] == 0
                SL[c] = hcat([V[:,j]; spzeros(n-nz, 1)], [spzeros(n-nz, 1); V[:,j]])
            else
                SL[c] = [spzeros(n-nz, 1); V[:,j]]
            end
        end
    end
    SU = copy(SL)
    for j = 1:length(SL)
        SL[j] = -1 / m / dx^2 * SL[j]
    end

    # nonlinear functions
    f = Vector(undef, length(brpts))
    for j = 1:interval-1
        coeff=brpts[j]
        f[j] = lambda -> exp(im * sqrt((m * (lambda - coeff * one(lambda)))))
    end
    for j = interval:length(brpts)
        coeff=brpts[j]
        f[j] = lambda -> exp(-sqrt((m * (-lambda + coeff * one(lambda)))))
    end

    # finally assemble nep instance; note that the nonlinear matrices are defined by their LU factors only
    C = map(k -> LowRankMatrixAndFunction(similar(SL[k], 0, 0), SL[k], SU[k], f[k]), 1:length(f))
    nep = SumNEP(PEP(B), LowRankFactorizedNEP(C))

    return nep, brpts, U0
end
