function particle_init(interval)
    nep, brpts, U0 = particle_nlep(interval)

    # interval
    sep = 1e-4
    if interval == 1
        xmin = -U0
        xmax = brpts[interval] - sep
        Xi = logspace(-6, 6, 10000) + brpts[interval]
    elseif interval > 1
        xmin = brpts[interval-1] + sep
        xmax = brpts[interval] - sep
        Xi1 = -logspace(-6, 6, 5000) + brpts[interval-1]
        Xi2 = logspace(-6, 6, 5000) + brpts[interval]
        Xi = [Xi1; Xi2]
    else
        error("Invalid interval: $interval")
    end

    # define target set Sigma
    Sigma = [xmin + 0im, xmax + 0im]

    # options
    srand(5)
    v0 = randn(nep.n)
    nodes = linspace(xmin, xmax, 11)
    nodes = nodes[2:2:end]
    funres = (Lam, X) -> particle_residual(Lam, X, nep)

    return nep, Sigma, Xi, v0, nodes, funres, xmin, xmax
end

function particle_nlep(interval)
    # MATLAB version author: William Vandenberghe

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

    x = kron(x_x, ones(z_z))
    z = kron(ones(x_x), z_z)

    w1 = 1 * nm
    w2 = 1.1 * nm
    l = 4 * nm
    U0 = 3 * eV

    # potential U
    U = zeros(x)
    U[abs.(z) .< w1] = -U0
    U[(abs.(z) .< w2) .& (abs.(x) .< l/2)] = -U0

    # m
    m = 0.2

    # branch points
    n = nx * nz
    Dxx_x = SymTridiagonal(ones(nx) * -2 / dx^2, ones(nx-1) / dx^2)
    Dzz_z = SymTridiagonal(ones(nz) * -2 / dx^2, ones(nz-1) / dz^2)
    H_L = -1/m*Dzz_z + diagm(U[1:nz])
    H_R = -1/m*Dzz_z + diagm(U[end-nz+(1:nz)])

    if norm(H_L - H_R) == 0
        # symmetric potential
        D,V = eig(H_L)
        i = sortperm(D)
        d = D[i]
        V = V[:,i]
        p = zeros(nz) # 0 (left+right)
    else
        # asymmetric potential
        D_L,V_L = eig(H_L)
        D_R,V_R = eig(H_R)
        V = [V_L; V_R]
        p = [-ones(nz); ones(nz)]
        DD = [D_L; D_R]
        i = sortperm(DD)
        d = DD[i]
        V = V[:,i]
        p = p[i] # -1 (left), 1 (right)
    end

    # matrix H
    H = -1/m * (kron(sparse(Dxx_x), speye(nz)) + kron(speye(nx), sparse(Dzz_z))) + spdiagm(U, 0, n, n)

    # branch points
    brpts = unique(d)

    # polynomial matrices
    B = [H, -speye(n)]

    # nonlinear matrices
    SL = Vector(length(brpts))
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
    S = Vector(length(brpts))
    for j = 1:length(S)
        S[j] = SL[j] * SU[j]'
    end

    # nonlinear functions
    f = Vector(length(brpts))
    for j = 1:interval-1
        f[j] = lambda -> expm(im * sqrtm(full(m * (lambda - brpts[j] * eye(lambda)))))
    end
    for j = interval:length(brpts)
        f[j] = lambda -> expm(-sqrtm(full(m * (-lambda + brpts[j] * eye(lambda)))))
    end

    # nlep
    NLEP = Dict("B" => B, "C" => S, "L" => SL, "U" => SU, "f" => f)

    # finally assemble nep instance
    C = map(k -> SPMFLowRankMatrix(S[k], SL[k], SU[k], f[k]), 1:length(S))
    nep = SPMFLowRankNEP(size(B[1], 1), B, C)

    return nep, brpts, U0
end

function particle_residual(Lambda, X, nep)
    function funA(lam, x)
        A = nep.B[1] * x
        for j = 2:length(nep.B)
            A += lam^(j-1) * (nep.B[j]*x)
        end
        as_matrix(x::Number) = (M = Matrix{eltype(x)}(1,1); M[1] = x; M)
        for j = 1:length(nep.C)
            A += nep.C[j].f(as_matrix(lam))[1] * (nep.C[j].A*x)
        end
        return A
    end

    map(i -> norm(funA(Lambda[i], X[:,i])), 1:length(Lambda))
end
