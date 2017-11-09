export jd_quad

using IterativeSolvers


jd_quad(nep::NEP;params...) = jd_quad(Complex128,nep;params...)
function jd_quad{T, T_orth<:IterativeSolvers.OrthogonalizationMethod}(::Type{T},
                    nep::ProjectableNEP;
                    orthmethod::Type{T_orth} = DGKS,
                    errmeasure::Function = default_errmeasure(nep::NEP),
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

    V::Array{T,2} = zeros(T, size(nep,1), maxit)
    V[:,1] = u
    dummy_vector = zeros(T,maxit)

    for k=1:maxit
        err=errmeasure(λ,u)
        if (displaylevel>0)
          println("Iteration:",k," errmeasure:",err)
        end
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

        u[:] = W*s

        Mdu::Array{T,1} = compute_Mlincomb(nep,λ,u,[1],1)
        P1 = eye(T,n) - Mdu*u'/(u'*Mdu)
        P2 = eye(T,n) - u*u'

        r = compute_Mlincomb(nep,λ,u)

        MP2 = zeros(T,size(nep,1),size(P2,2))
        for ii = 1:size(P2,2)
            MP2[:,ii] = compute_Mlincomb(nep,λ,P2[:,ii],[1],0)
        end
        X = P1*MP2

        # Least squares, -pseudo_inv(X)*r
        Q,R = qr(X)
        w[:] = -R\(Q'*r)

        # orthogonalization
        orthogonalize_and_normalize!(W, w, dummy_vector[1:k], orthmethod)

        println("Iteration: ",k," norm of residual:", compute_resnorm(nep,λ,u))
    end

    err=errmeasure(λ,u)
    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(λ,v,err,msg))
end


function jd_inner_eig_solver(::Type{T_orig_nep}, T, proj_nep, eigsolvertype) where {T_orig_nep <: ProjectableNEP}
    error("NOT IMPLEMTED YET") #TODO: Implement this!
end

function jd_inner_eig_solver(::Type{T_orig_nep}, T, proj_nep, eigsolvertype) where {T_orig_nep <: PEP}
    pep_temp = PEP(get_Av(proj_nep))
    Dc,Vc = polyeig(T,pep_temp,eigsolvertype)
    c = sortperm(abs.(Dc))
    λ = Dc[c[1]]
    s = Vc[:,c[1]]
    return λ, s
end
