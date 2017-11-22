export jd

using IterativeSolvers


jd(nep::NEP;params...) = jd(Complex128,nep;params...)
function jd{T, T_orth<:IterativeSolvers.OrthogonalizationMethod}(::Type{T},
                    nep::ProjectableNEP;
                    orthmethod::Type{T_orth} = DGKS,
                    errmeasure::Function = default_errmeasure(nep::NEP),
                    linsolvercreator::Function=default_linsolvercreator,
                    proj_eig_solver::Function = default_proj_eig_solver(nep::ProjectableNEP),
                    Neig = 1,
                    tol = eps(real(T))*100,
                    maxit = 100,
                    λ = zero(T),
                    v0 = randn(size(nep,1)),
                    displaylevel = 0)

    n = size(nep,1)
    if (maxit > n)
        warn("maxit = ", maxit, " is larger than size of NEP = ", n,". Setting maxit = size(nep,1)")
        maxit = n
    end

    λ_vec::Array{T,1} = Array{T,1}(Neig)
    u_vec::Array{T,2} = Array{T,2}(n,Neig)
    u::Array{T,1} = Array{T,1}(v0)/norm(v0)

    proj_nep = create_proj_NEP(nep)

    V::Array{T,2} = zeros(T, size(nep,1), maxit+1)
    V[:,1] = u
    dummy_vector = zeros(T,maxit+1)
    conveig = 0

    for k=1:maxit
        err=errmeasure(λ,u)
        @ifd(print("Iteration: ", k, " converged eigenvalues: ", conveig, " errmeasure: ", err, "\n"))
        if convergence_criterion(err, tol, λ, λ_vec, conveig, T)
            conveig += 1
            λ_vec[conveig] = λ
            u_vec[:,conveig] = u
            if (conveig == Neig)
                return (λ_vec,u_vec)
            end
        end

        # Projected matrices
        W = view(V, :, 1:k) # extact subarrays, memory-CPU efficient
        w = view(V, :, k+1) # next vector position
        set_projectmatrices!(proj_nep, W, W)

        # find the eigenvalue with smallest absolute value of projected NEP
        λ,s = proj_eig_solver(typeof(nep), T, proj_nep, conveig+1)
        s = s/norm(s)

        # the approximate eigenvector
        u[:] = W*s

        # solve for the basis extension using equation (9) to avoid matrix access.
        # the vector should be orthogonal to u
        pk::Array{T,1} = compute_Mlincomb(nep,λ,u,[1],1)
        linsolver = linsolvercreator(nep,λ)
        tt::Array{T,1} = lin_solve(linsolver, pk, tol=tol) # M(λ)\pk
        α = norm(u)^2/dot(tt,u)
        w[:] = -u + α*tt

        # orthogonalization
        orthogonalize_and_normalize!(W, w, dummy_vector[1:k], orthmethod)
    end

    err=errmeasure(λ,u)
    # Check one last time since space was expanded since prevous
    if convergence_criterion(err, tol, λ, λ_vec, conveig, T)
        conveig += 1
        λ_vec[conveig] = λ
        u_vec[:,conveig] = u
        if (conveig == Neig)
            return (λ_vec,u_vec)
        end
    end

    msg="Number of iterations exceeded. maxit=$(maxit) and only $(conveig) eigenvalues converged out of $(Neig)."
    throw(NoConvergenceException(cat(1,λ_vec[1:conveig],λ),cat(2,u_vec[:,1:conveig],u),err,msg))
end


function convergence_criterion(err, tol, λ, λ_vec, conveig, T)
    #Small error and not already found (expception if it is the first)
    return (err < tol) && (conveig == 0 || all( abs.(λ - λ_vec[1:conveig])./abs.(λ_vec[1:conveig]) .> eps(real(T))*1e5 ) )
end

function default_proj_eig_solver(nep::ProjectableNEP)
    f=function(::Type{T_orig_nep}, T, proj_nep, N) where {T_orig_nep <: ProjectableNEP}
        defaul_jd_inner_eig_solver(T_orig_nep, T, proj_nep, N)
    end
    return f
end

function defaul_jd_inner_eig_solver(::Type{T_orig_nep}, T, proj_nep, N) where {T_orig_nep <: ProjectableNEP}
    λ,s = sgiter(T, proj_nep, N)
    return λ, s
end

function defaul_jd_inner_eig_solver(::Type{T_orig_nep}, T, proj_nep, N) where {T_orig_nep <: PEP}
    pep_temp = PEP(get_Av(proj_nep))
    Dc,Vc = polyeig(T,pep_temp)
    c = sortperm(abs.(Dc))
    λ = Dc[c[N]]
    s = Vc[:,c[N]]
    return λ, s
end
