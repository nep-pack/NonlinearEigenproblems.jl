



export jd_quad

jd_quad(nep::NEP;params...) = jd_quad(Complex128,nep;params...)
function jd_quad{T}(::Type{T},
                    nep::ProjectableNEP;
                    errmeasure::Function =
                    default_errmeasure(nep::NEP),
                    tolerance=eps(real(T))*100,
                    maxit=100,
                    λ=zero(T),
                    v0=randn(size(nep,1)),
                    displaylevel=0,
                    eigsolvertype::DataType=DefaultEigSolver)


    λ::T = T(λ)
    v0 = Array{T,1}(v0)
    v::Array{T,1} = v0/norm(v0)
    n = size(nep,1)

    proj_nep = create_proj_NEP(nep)

    V::Array{T,2} = zeros(T,size(nep,1),1)
    V[:,1] = v
    u::Array{T,1} = v

    #loop...
    for kk=1:maxit

        err=errmeasure(λ,u)
        if (displaylevel>0)
          println("Iteration:",k," errmeasure:",err)
        end
        if (err< tolerance)
            return (λ,u)
        end

        # Projected matrices
        set_projectmatrices!(proj_nep,V,V)

        # Create the projected NEP problem and
        # find the eigenvalue with smallest absolute value
        Dc,Vc = jd_inner_eig_solver(typeof(nep), T, proj_nep, eigsolvertype)
        c = sortperm(abs.(Dc))

        λ = Dc[c[1]]
        s = Vc[:,c[1]]
        s = s/norm(s)

        u = V*s

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
        t = -R\(Q'*r)

        #Modified Gram-Schmidt
        for ii=1:kk
            temp = dot(V[:,ii],t)
            t += - temp*V[:,ii]
        end
        # reorthogonalization
        for ii=1:kk
              temp = dot(V[:,ii],t)
              t += - temp*V[:,ii]
        end
        v = t/norm(t);

        # Update the search space V
        V = [V v]

        println("Iteration: ",kk," norm of residual:", compute_resnorm(nep,λ,u))
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
    return Dc,Vc
end
