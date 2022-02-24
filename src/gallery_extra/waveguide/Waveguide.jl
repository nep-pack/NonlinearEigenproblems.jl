
###########################################################
# Waveguide eigenvalue problem (WEP)
# Sum of products of matrices and functions (SPMF)
"""
Waveguide eigenvalue problem (WEP)
Sum of products of matrices and functions (SPMF)
"""
function assemble_waveguide_spmf_fd(nx::Integer, nz::Integer, hx, Dxx::SparseMatrixCSC, Dzz::SparseMatrixCSC, Dz::SparseMatrixCSC, C1::SparseMatrixCSC, C2T::SparseMatrixCSC, K::Union{Matrix{ComplexF64},Matrix{Float64}}, Km, Kp, pre_Schur_fact::Bool)
    Ix = sparse(ComplexF64(1)I, nx, nx)
    Iz = sparse(ComplexF64(1)I, nz, nz)
    Q0 = kron(Ix, Dzz) + kron(Dxx, Iz) + sparse(Diagonal(vec(K)))
    Q1 = kron(Ix, 2*Dz)
    Q2 = kron(Ix, Iz)

    A = Vector{SparseMatrixCSC}(undef, 3+2*nz)
    A[1] = hvcat((2,2), Q0, C1, C2T, spzeros(ComplexF64,2*nz, 2*nz) )
    A[2] = hvcat((2,2), Q1, spzeros(ComplexF64,nx*nz, 2*nz), spzeros(ComplexF64,2*nz, nx*nz), spzeros(ComplexF64,2*nz, 2*nz) )
    A[3] = hvcat((2,2), Q2, spzeros(ComplexF64,nx*nz, 2*nz), spzeros(ComplexF64,2*nz, nx*nz), spzeros(ComplexF64,2*nz, 2*nz) )

    f = Vector{Function}(undef, 3+2*nz)
    f[1] = λ -> one(λ)
    f[2] = λ -> λ
    f[3] = λ -> λ^2

    R, Rinv = generate_R_matvecs(nz)
    S = generate_S_function(nz, hx, Km, Kp)

    for j = 1:nz
        f[j+3] = λ -> S(λ, j)
        e_j = zeros(nz)
        e_j[j] = 1
        E_j = [R(e_j); spzeros(ComplexF64,nz)]
        E_j = E_j * (E_j/nz)'
        A[j+3] =  hvcat((2,2), spzeros(ComplexF64,nx*nz,nx*nz), spzeros(ComplexF64,nx*nz, 2*nz), spzeros(ComplexF64,2*nz, nx*nz), E_j)
    end
    for j = 1:nz
        f[j+nz+3] = λ -> S(λ, nz+j)
        e_j = zeros(nz)
        e_j[j] = 1
        E_j = [spzeros(ComplexF64,nz); R(e_j)]
        E_j = E_j * (E_j/nz)'
        A[j+nz+3] =  hvcat((2,2), spzeros(ComplexF64,nx*nz,nx*nz), spzeros(ComplexF64,nx*nz, 2*nz), spzeros(ComplexF64,2*nz, nx*nz), E_j)
    end
    return SPMF_NEP(A,f;Schur_fact=pre_Schur_fact)
end


###########################################################
# Generate R-matrix
# Part of defining the P-matrix, see above, Jarlebring-(1.6) and Ringh-(2.8) and Remark 1
# OBS: R' = nz * Rinv, as noted in Ringh between (2.8) and Remark 1
function generate_R_matvecs(nz)
    # The scaled FFT-matrix R
    p = (nz-1)/2
    bb = exp.(-2im*pi*((1:nz).-1)*(-p)/nz)  # scaling to do after FFT
    function R(X) # Note! Only works for vectors or one-dim matrices
        return reverse(bb .* fft(vec(X)), dims = 1)
    end
    bbinv = 1 ./ bb # scaling to do before inverse FFT
    function Rinv(X)
        return ifft(bbinv .* reverse(vec(X), dims = 1))
    end
    return R, Rinv
end


###########################################################
# Generate S-function for matrix argument
# Part of defining the P-matrix, see above, Jarlebring-(2.4) and Ringh-(2.8)(2.3)
function generate_S_function(nz, hx, Km, Kp)
    # Constants from the problem
    p = (nz-1)/2
    d0 = -3/(2*hx)
    b = 4*pi*1im * (-p:p)
    cM = Km^2 .- 4*pi^2 * ((-p:p).^2)
    cP = Kp^2 .- 4*pi^2 * ((-p:p).^2)

    # Note: λ should be scalar or matrix (not vector)
    betaM = function(λ, j)
        return λ^2 + b[j]*λ + cM[j]*one(λ)
    end
    betaP = function(λ, j)
        return λ^2 + b[j]*λ + cP[j]*one(λ)
    end

    sM = function(λ, j)
        return  1im*sqrt_schur_pos_imag(betaM(λ, j)) + d0*one(λ)
    end
    sP = function(λ, j)
        return  1im*sqrt_schur_pos_imag(betaP(λ, j)) + d0*one(λ)
    end


    S = function(λ, j)
        if j <= nz
            return sM(λ,j)
        elseif j <= 2*nz
            return sP(λ,(j-nz))
        else
            error("The chosen j = ", j, "but the setup nz = ", nz, ". Hence j>2nz which is illegal.")
        end
    end

    return S
end


###########################################################
"""
    sqrt_schur_pos_imag(A::AbstractMatrix)
 Computes the matrix square root on the 'correct branch',
 that is, with positivt imaginary part. Similar to Schur method
 in Algorithm 6.3 in Higham matrix functions.
"""
sqrt_schur_pos_imag(A::Number) = sqrt_pos_imag(A)
function sqrt_schur_pos_imag(A::AbstractMatrix)
    n = size(A,1);
    AA = Matrix{ComplexF64}(A);
    (T, Q, ) = schur(AA)
    U = zeros(ComplexF64,n,n);
    for i = 1:n
        U[i,i] = sqrt_pos_imag(T[i,i])
    end
    private_inner_loops_sqrt!(n, U, T)
    return Q*U*Q'
end
#Helper function executing the inner loop (more Juliaesque)
function private_inner_loops_sqrt!(n, U, T)
    temp = zero(ComplexF64);
    for j = 2:n
        for i = (j-1):-1:1
            temp *= zero(ComplexF64);
            for k = (i+1):(j-1)
                temp += U[i,k]*U[k,j]
            end
            U[i,j] = (T[i,j] - temp)/(U[i,i]+U[j,j])
        end
    end
end

"""
    sqrt_pos_imag(a::ComplexF64) and sqrt_pos_imag(a::Float64)
 Helper function: Computes the scalar square root on the 'correct branch',
 that is, with positivt imaginary part.
"""
function sqrt_pos_imag(a::ComplexF64)
    imag_sign = sign(imag(a))
    if imag_sign == 0 #Real value in complex type
        sqrt(a)
    else
        sign(imag(a))*sqrt(a)
    end
end
function sqrt_pos_imag(a::Float64)
    return sqrt(a)
end

    function Pinv(nep,λ,x)
        return vec(  [R(nep,Rinv(nep,x[1:Int(end/2)]) ./ sM(nep,λ));
                      R(nep,Rinv(nep,x[Int(end/2)+1:end]) ./ sP(nep,λ))  ]  )
    end

    function R(nep,X) # Note! Only works for vectors or one-dim matrices
        return reverse(nep.bb .* fft(vec(X)), dims = 1)
    end

    function Rinv(nep,X)
        return ifft(nep.bbinv .* reverse(vec(X), dims = 1))
    end

    function betaM(nep,λ)
        return λ^2 .+ nep.b*λ + nep.cM
    end

    function betaP(nep,λ)
        return λ^2 .+ nep.b*λ + nep.cP
    end

    function sM(nep,λ)
        bbeta = betaM(nep,λ)
        return 1im*sign.(imag(bbeta)).*sqrt.(bbeta) .+ nep.d0
    end

    function sP(nep,λ)
        bbeta = betaP(nep,λ)
        return 1im*sign.(imag(bbeta)).*sqrt.(bbeta) .+ nep.d0
    end

###########################################################
# Waveguide eigenvalue problem - WEP
# A more optimized (native) implementation of the WEP with FD discretization
"""
    WEP_FD

Waveguide eigenvalue problem

A more optimized implementation of the WEP for FD-discretization.\\
Closer to what is done in the article:
''E. Ringh, and G. Mele, and J. Karlsson, and E. Jarlebring,
Sylvester-based preconditioning for the waveguide eigenvalue problem,
Linear Algebra and its Applications''
"""
    struct WEP_FD <: WEP
        nx::Int
        nz::Int
        hx::Float64
        hz::Float64
        Dxx::SparseMatrixCSC{Float64,Int}
        Dzz::SparseMatrixCSC{Float64,Int}
        Dz::SparseMatrixCSC{Float64,Int}
        Iz::SparseMatrixCSC{Float64,Int}
        C1::SparseMatrixCSC{Float64,Int}
        C2T::SparseMatrixCSC{Float64,Int}
        k_bar::ComplexF64
        K::Matrix{ComplexF64}
        p::Integer
        d0::Float64
        d1::Float64
        d2::Float64
        b::Vector{ComplexF64}
        cM::Vector{ComplexF64}
        cP::Vector{ComplexF64}
        bb::Vector{ComplexF64}
        bbinv::Vector{ComplexF64}

        function WEP_FD(nx, nz, hx, hz, Dxx, Dzz, Dz, C1, C2T, K, Km, Kp)
            n = nx*nz + 2*nz
            k_bar = mean(K)
            K_scaled = K-k_bar*ones(ComplexF64,nz,nx)

            p = (nz-1)/2

            d0 = -3/(2*hx)
            d1 = 2/hx
            d2 = -1/(2*hx)

            b = 4*pi*1im * (-p:p)
            cM = Km^2 .- 4*pi^2 * ((-p:p).^2)
            cP = Kp^2 .- 4*pi^2 * ((-p:p).^2)

            bb = exp.(-2im*pi*((1:nz).-1)*(-p)/nz)  # scaling to do after FFT
            bbinv = 1 ./ bb # scaling to do before inverse FFT

            Iz = sparse(ComplexF64(1)I, nz, nz)

            return new(nx, nz, hx, hz, Dxx, Dzz, Dz, Iz, C1, C2T, k_bar, K_scaled, p, d0, d1, d2, b, cM, cP, bb, bbinv)
        end
    end

    function A(nep::WEP_FD, λ, d=0)
        nz = nep.nz
        if(d == 0)
            return nep.Dzz + 2*λ*nep.Dz + (λ^2 + nep.k_bar)*nep.Iz
        elseif(d == 1)
            return 2*nep.Dz + 2*λ*nep.Iz
        elseif(d == 2)
            return 2*nep.Iz
        else
            return spzeros(ComplexF64, nz, nz)
        end
    end

    function B(nep::WEP_FD, λ, d=0)
        if(d == 0)
            return nep.Dxx
        else
            return spzeros(ComplexF64, nx, nx)
        end
    end

    #Inverses of the boundary operators, Ringh - (2.8)
    #To be used in the Schur-complement- and SMW-context.
    function P_inv_m(nep::WEP_FD, λ, v)
        nz = nep.nz
        coeffs = zeros(ComplexF64, nz)
            aa = 1.0
        for j = 1:nz
            bb = nep.b[j]
            cc = nep.cM[j]
            coeffs[j] = 1im*sqrt_derivative(aa, bb, cc, 0, λ) + nep.d0
        end
        return R(nep,Rinv(nep,v) ./ coeffs)
    end

    function P_inv_p(nep::WEP_FD, λ, v)
        nz = nep.nz
        coeffs = zeros(ComplexF64, nz)
            aa = 1.0
        for j = 1:nz
            bb = nep.b[j]
            cc = nep.cP[j]
            coeffs[j] = 1im*sqrt_derivative(aa, bb, cc, 0, λ) + nep.d0
        end
        return R(nep,Rinv(nep,v) ./ coeffs)
    end


    function size(nep::WEP_FD)
        n = nep.nx*nep.nz + 2*nep.nz
        return (n,n)
    end
    function size(nep::WEP_FD, dim)
        n = nep.nx*nep.nz + 2*nep.nz
        return n
    end


    function issparse(nep::WEP_FD)
        return false
    end





"""
    compute_Mlincomb(nep::WEP_FD, λ::Number, V; a=ones(ComplexF64,size(V,2)))

Specialized for Waveguide Eigenvalue Problem discretized with Finite Difference\\\\
 Computes the linear combination of derivatives\\
 ``Σ_i a_i M^{(i)}(λ) v_i``
"""
    function compute_Mlincomb(nep::WEP_FD, λ::Number, V::AbstractVecOrMat,
                              a::Vector=ones(ComplexF64,size(V,2)))
        na = size(a,1)
        nv = size(V,1)
        mv = size(V,2)
        n_nep = size(nep,1)
        if(na != mv)
            error("Incompatible sizes: Number of coefficients = ", na, ", number of vectors = ", mv, ".")
        end
        if(nv != n_nep)
            error("Incompatible sizes: Length of vectors = ", nv, ", size of NEP = ", n_nep, ".")
        end
        nx = nep.nx
        nz = nep.nz
        max_d = na - 1 #Start on 0:th derivative

        V1 = view(V, 1:nx*nz, :)
        V1_mat = reshape(V1, nz, nx, na)
        V2 = view(V, nx*nz+1:n_nep, :)

        # Compute the top part (nx*nz)
        y1_mat::Matrix{ComplexF64} = (A(nep,λ) * V1_mat[:,:,1] + V1_mat[:,:,1] * B(nep,λ)  +  nep.K .* V1_mat[:,:,1])*a[1]
        for d = 1:min(max_d,3)
            y1_mat[:,:] += A(nep,λ,d) * V1_mat[:,:,d+1] * a[d+1];
        end
        y1::Vector{ComplexF64} = y1_mat[:]
        y1[:] += nep.C1 * V2[:,1] * a[1]

        # Compute the bottom part (2*nz)
        D::Matrix{ComplexF64} = zeros(ComplexF64, 2*nz, na)
        cMP = vcat(nep.cM, nep.cP)
        for j = 1:2*nz
            aa = 1
            bb = nep.b[rem(j-1,nz)+1]
            cc = cMP[j]
            der_coeff = 1im*sqrt_derivative(aa, bb, cc, max_d, λ)
            for jj = 1:na
                D[j, jj] = der_coeff[jj]
            end
        end

        #Multiplication with diagonal matrix optimized by working "elementwise" Jarlebring-(4.6)
        y2_temp::Vector{ComplexF64} =
            (D[:,1] .+ nep.d0) .* [Rinv(nep,V2[1:nz,1]);
                                   Rinv(nep,V2[nz+1:2*nz,1])]*a[1]
        for jj = 2:na
            #Multiplication with diagonal matrix optimized by working "elementwise" Jarlebring-(4.6)
            y2_temp[:] += D[:,jj] .* [Rinv(nep,V2[1:nz,jj]);
                                      Rinv(nep,V2[nz+1:2*nz,jj])] *a[jj]
        end
        y2::Vector{ComplexF64} = [R(nep,y2_temp[1:nz,1]);
                                  R(nep,y2_temp[nz+1:2*nz,1])]
        y2[:] += nep.C2T * V1[:,1]*a[1] #Action of C2T. OBS: Add last because of implcit storage in R*D_i*R^{-1}*v_i

        return vcat(y1, y2)
    end

"""
    compute_Mder(nep::WEP,λ::Number,i::Integer=0)
"""
    function compute_Mder(nep::WEP,λ::Number,i::Integer=0)
        error("The WEP does not implement this function. If this was called in a situation where you want to solve linear systems please look at `WEPLinSolverCreator`")
    end


###########################################################
# Linear Solvers for WEP

    # Matrix vector operations for the Schur complement (to be used in GMRES call)
    # Matrix-vector product according to Ringh (2.13) and (3.3)
    struct SchurMatVec
        nep::WEP_FD
        λ::ComplexF64
        SchurMatVec(nep::WEP_FD, λ::Union{ComplexF64,Float64}) = new(nep, λ)
    end

    function *(M::SchurMatVec,v::AbstractVector)
        λ = M.λ
        nep = M.nep

        X = reshape(v, nep.nz, nep.nx)
        return vec(  vec( A(nep,λ)*X + X*B(nep,λ) + nep.K.*X ) - nep.C1 * Pinv(nep, λ, nep.C2T*v)  )
    end

    function (M::SchurMatVec)(v::AbstractVector) #Overload the ()-function so that a SchurMatVec struct can act and behave like a function
        return M*v
    end

    function size(M::SchurMatVec, dim=-1)
        n = M.nep.nx*M.nep.nz
        if (dim==-1)
            return (n,n)
        else
            return n
        end
    end

    function eltype(M::SchurMatVec)
        return ComplexF64
    end


    # GMRES Solver
    struct WEPGMRESLinSolver<:LinSolver
        schur_comp::LinearMap{ComplexF64}
        kwargs
        gmres_log::Bool
        nep::WEP_FD
        λ::ComplexF64

        function WEPGMRESLinSolver(nep::WEP_FD,λ::Union{ComplexF64,Float64},kwargs)
            f = SchurMatVec(nep, λ)
            schur_comp = LinearMap{ComplexF64}(f, nep.nx*nep.nz, ismutating=false, issymmetric=false, ishermitian=false);
            gmres_log = false
            for elem in kwargs
                gmres_log |= ((elem[1] == :log) && elem[2])
            end
            return new(schur_comp, kwargs, gmres_log,nep,λ)
        end
    end

    function WEP_inner_lin_solve(solver::WEPGMRESLinSolver, rhs::Vector, reltol)
        if( solver.gmres_log )
            q, convhist = gmres(solver.schur_comp, rhs; reltol=reltol, solver.kwargs...)
        else
            q = gmres(solver.schur_comp, rhs; reltol=reltol, solver.kwargs...)
        end
        return q
    end


    # Direct Backslash solver
    struct WEPBackslashLinSolver<:LinSolver
        schur_comp::SparseMatrixCSC{ComplexF64,Int}
        nep::WEP_FD
        λ::ComplexF64

        function WEPBackslashLinSolver(nep::WEP_FD, λ::Union{ComplexF64,Float64}, kwargs=())
            schur_comp = construct_WEP_schur_complement(nep, λ)
            return new(schur_comp, nep, λ)
        end
    end

    function WEP_inner_lin_solve(solver::WEPBackslashLinSolver, rhs::Vector, tol)
        return solver.schur_comp \ rhs
    end


    # Direct pre-factorized solver
    struct WEPFactorizedLinSolver<:LinSolver
        schur_comp_fact
        nep::WEP_FD
        λ::ComplexF64

        function WEPFactorizedLinSolver(nep::WEP_FD, λ::Union{ComplexF64,Float64}, kwargs=())
            schur_comp_fact = factorize(construct_WEP_schur_complement(nep, λ))
            return new(schur_comp_fact, nep, λ)
        end
    end

    function WEP_inner_lin_solve(solver::WEPFactorizedLinSolver, rhs::Vector, tol)
        return solver.schur_comp_fact \ rhs
    end

"""
    struct WEPLinSolverCreator <: LinSolverCreator
    function WEPLinSolverCreator(;solver_type=:backslash,kwargs=())

The a linsolver creator for the waveguide eigenvalue problem. The
kwarg `solver_type` is one of `:backslash`, `:factorized`, `:gmres`.
The `kwargs` keyword argument is passed to the solver.
"""
    struct WEPLinSolverCreator <: LinSolverCreator
        solver_type::Symbol
        kwargs;
        function WEPLinSolverCreator(;solver_type=:factorized,kwargs=())
            return new(solver_type,kwargs);
        end
    end

    function create_linsolver(creator::WEPLinSolverCreator,nep,λ)
        if (!(nep isa WEP_FD))
            t=typeof(nep)
            error("WEPLinSolver can only be used in combination with WEPs: type(nep)=$t");
        end
        if (creator.solver_type==:backslash)
            return WEPBackslashLinSolver(nep, λ, creator.kwargs)
        elseif (creator.solver_type==:gmres)
            return WEPGMRESLinSolver(nep, λ, creator.kwargs)
        elseif (creator.solver_type==:factorized)
            return WEPFactorizedLinSolver(nep, λ, creator.kwargs)
        else
            s=creator.solver_type
            error("Unknown type of solver_type in linsolvercreator:$s");
        end
    end

    # Helper functions for WEP LinSolvers. To avoid code repetition.
    # Assembls the full Schur-complement, used in both Backslash and LU solvers
    function construct_WEP_schur_complement(nep::WEP_FD, λ::Union{ComplexF64,Float64})
        nz = nep.nz
        nx = nep.nx
        Inz = sparse(ComplexF64(1)I, nz, nz)
        Inx = sparse(ComplexF64(1)I, nx, nx)

        Pinv_minus = Matrix{ComplexF64}(undef, nz, nz)
        Pinv_plus = Matrix{ComplexF64}(undef, nz, nz)
        e = zeros(ComplexF64,nz)
        for i = 1:nz
            e[:] .= 0
            e[i] = 1
            Pinv_minus[:,i] = P_inv_m(nep, λ, e)
            Pinv_plus[:,i] = P_inv_p(nep, λ, e)
        end

        E = spzeros(nx,nx)
        E[1,1] = nep.d1/(nep.hx^2)
        E[1,2] = nep.d2/(nep.hx^2)
        EE = spzeros(nx,nx)
        EE[nx,nx] = nep.d1/(nep.hx^2)
        EE[nx,nx-1] = nep.d2/(nep.hx^2)

        # Kronecker product form of Ringh - Proposition 3.1
        return kron(copy(B(nep,λ)'), Inz) + kron(Inx, A(nep,λ)) + sparse(Diagonal(nep.K[:])) - kron(E, Pinv_minus) - kron(EE, Pinv_plus)
    end

    # lin_solve function to wrapp all the WEP linear solvers.
    # Since Schur-complement transformations are the same.
    # Does transforming between that and the full system.
    # Ringh - Proposition 2.1, see also Algorithm 2, step 10-11.
    function lin_solve(solver::Union{WEPBackslashLinSolver,WEPGMRESLinSolver,WEPFactorizedLinSolver}, x::Vector; tol=eps(Float64))
    # Ringh - Proposition 2.1
        λ = solver.λ
        nep = solver.nep

        x_int = x[1:(nep.nx*nep.nz)]
        x_ext = x[((nep.nx*nep.nz)+1):((nep.nx*nep.nz) + 2*nep.nz)]
        rhs =  vec(  x_int - nep.C1*Pinv(nep, λ, x_ext))

        q = WEP_inner_lin_solve(solver, rhs, tol)

        return [q; vec(Pinv(nep, λ, -nep.C2T * q + x_ext))]

    end


###########################################################
#Square root of second degree polynomial (Gegenbauer polynomials)
#Jarlebring - Appendix C
"""
    sqrt_derivative(a,b,c, d=0, x=0)

Computes all d derivatives of sqrt(a*z^2 + b*z + c)
in the point z = x.
Returns a d+1 vector with all numerical values
"""
function sqrt_derivative(a,b,c, d=0, x=0)
    if(d<0)
        error("Cannot take negative derivative. d = ", d)
    end

    aa = a
    bb = b + 2*a*x
    cc = c + a*x^2 + b*x

    derivatives = zeros(ComplexF64,d+1)

    yi = sqrt_pos_imag(cc)
    derivatives[1] = yi
    if( d==0 )
        return derivatives[1] #OBS: If only function is sought, return it as an integer and not array
    end

    yip1 = bb/(2*sqrt_pos_imag(cc))
    fact = Float64(1)
    derivatives[2] = yip1 * fact
    if( d==1 )
        return derivatives
    end

    yip2 = zero(ComplexF64)
    for i = 2:d
        m = i - 2
        yip2 = - (2*aa*(m-1)*yi  +  bb*(1+2*m)*yip1) / (2*cc*(2+m))
        fact *= i

        yi = yip1
        yip1 = yip2

        derivatives[i+1] = yip2 * fact
    end
    return derivatives
end
