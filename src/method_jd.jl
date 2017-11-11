export jd

using IterativeSolvers


jd(nep::NEP;params...) = jd(Complex128,nep;params...)
function jd{T, T_orth<:IterativeSolvers.OrthogonalizationMethod}(::Type{T},
                    nep::ProjectableNEP;
                    orthmethod::Type{T_orth} = DGKS,
                    errmeasure::Function = default_errmeasure(nep::NEP),
                    linsolvercreator::Function=default_linsolvercreator,
                    tol = eps(real(T))*100,
                    maxit = 100,
                    λ = zero(T),
                    v0 = randn(size(nep,1)),
                    displaylevel = 0,
                    eigsolvertype::DataType = DefaultEigSolver)

    λ::T = T(λ)
    u::Array{T,1} = Array{T,1}(v0)/norm(v0)
    n = size(nep,1)

    proj_nep = create_proj_NEP(nep)

    V::Array{T,2} = zeros(T, size(nep,1), maxit+1)
    V[:,1] = u
    dummy_vector = zeros(T,maxit+1)

    for k=1:maxit
        err=errmeasure(λ,u)
        @ifd(print("Iteration:", k, " errmeasure:", err, "\n"))
        if (err< tol)
            return (λ,u)
        end

        # Projected matrices
        W = view(V, :, 1:k) # extact subarrays, memory-CPU efficient
        w = view(V, :, k+1); # next vector position
        set_projectmatrices!(proj_nep, W, W)

        # find the eigenvalue with smallest absolute value of projected NEP
        λ,s = jd_inner_eig_solver(typeof(nep), T, proj_nep, eigsolvertype)
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
    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(λ,u,err,msg))
end


function jd_inner_eig_solver(::Type{T_orig_nep}, T, proj_nep, eigsolvertype) where {T_orig_nep <: ProjectableNEP}
    λ,s = sgiter(T, proj_nep, 1)
    return λ, s
end

function jd_inner_eig_solver(::Type{T_orig_nep}, T, proj_nep, eigsolvertype) where {T_orig_nep <: PEP}
    pep_temp = PEP(get_Av(proj_nep))
    Dc,Vc = polyeig(T,pep_temp,eigsolvertype)
    c = sortperm(abs.(Dc))
    λ = Dc[c[1]]
    s = Vc[:,c[1]]
    return λ, s
end
