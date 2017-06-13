module Waveguide

using NEPCore
using NEPTypes
using LinSolvers
using IterativeSolvers


export gallery_waveguide
export generate_preconditioner



# Specializalized NEPs
export WEP_FD

# Preconditioner
export generate_smw_matrix


# We overload these
import NEPCore.compute_Mder
export compute_Mder
import NEPCore.compute_Mlincomb
export compute_Mlincomb
import LinSolvers.Mlincomb_matvec
import LinSolvers.lin_solve
export lin_solve

import Base.size
export size
import Base.issparse
export issparse
import Base.*
export *


include("waveguide_debug.jl")
include("waveguide_FD.jl")
include("waveguide_FEM.jl")
# OBS: include("waveguide_preconditioner.jl") is at the bottom so that needed types are declared before the include



"
Creates the NEP associated with example in

E. Ringh, and G. Mele, and J. Karlsson, and E. Jarlebring, 
Sylvester-based preconditioning for the waveguide eigenvalue problem,
Linear Algebra and its Applications

and

E. Jarlebring, and G. Mele, and O. Runborg
The waveguide eigenvalue problem and the tensor infinite Arnoldi method
SIAM J. Sci. Comput., 2017
"
function gallery_waveguide( nx::Integer = 3*5*7, nz::Integer = 3*5*7, waveguide::String = "TAUSCH", discretization::String = "FD", NEP_format::String = "SPMF",  delta::Number = 0.1)
    waveguide = uppercase(waveguide)
    NEP_format = uppercase(NEP_format)
    discretization = uppercase(discretization)
    if !isodd(nz)
        error("Variable nz must be odd! You have used nz = ", nz, ".")
    end


    # Generate the matrices and formulate the problem is the sought format
    if (NEP_format == "SPMF") && (discretization == "FD")
        K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        nep = assemble_waveguide_spmf_fd(nx, nz, hx, Dxx, Dzz, Dz, C1, C2T, K, Km, Kp)

    elseif (NEP_format == "WEP") && (discretization == "FD")
        K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        nep = WEP_FD(nx, nz, hx, hz, Dxx, Dzz, Dz, C1, C2T, K, Km, Kp)

    else
        error("The NEP-format '", NEP_format, "' is not supported for the discretization '", discretization, "'.")
    end

    println("Waveguide generated")
    return nep
end 


###########################################################
# Waveguide eigenvalue problem (WEP)
# Sum of products of matrices and functions (SPMF)
"""
 Waveguide eigenvalue problem (WEP)
Sum of products of matrices and functions (SPMF)
"""
function assemble_waveguide_spmf_fd(nx::Integer, nz::Integer, hx, Dxx::SparseMatrixCSC, Dzz::SparseMatrixCSC, Dz::SparseMatrixCSC, C1::SparseMatrixCSC, C2T::SparseMatrixCSC, K::Union{Array{Complex128,2},Array{Float64,2}}, Km, Kp)
    Ix = speye(nx,nx)
    Iz = speye(nz,nz)
    Q0 = kron(Ix, Dzz) + kron(Dxx, Iz) + spdiagm(vec(K))
    Q1 = kron(Ix, 2*Dz)
    Q2 = kron(Ix, Iz)

    A = Array(SparseMatrixCSC,3+2*nz)
    A[1] = hvcat((2,2), Q0, C1, C2T, spzeros(2*nz, 2*nz) )
    A[2] = hvcat((2,2), Q1, spzeros(nx*nz, 2*nz), spzeros(2*nz, nx*nz), spzeros(2*nz, 2*nz) )
    A[3] = hvcat((2,2), Q2, spzeros(nx*nz, 2*nz), spzeros(2*nz, nx*nz), spzeros(2*nz, 2*nz) )

    f = Array(Function, 3+2*nz)
    f[1] = λ -> eye(λ)
    f[2] = λ -> λ
    f[3] = λ -> λ^2

    R, Rinv = generate_R_matrix(nz)
    S = generate_S_function(nz, hx, Km, Kp)

    for j = 1:nz
        f[j+3] = λ -> S(λ, j)
        e_j = zeros(nz)
        e_j[j] = 1
        E_j = [R(e_j); spzeros(nz)]
        E_j = E_j * (E_j/nz)'
        A[j+3] =  hvcat((2,2), spzeros(nx*nz,nx*nz), spzeros(nx*nz, 2*nz), spzeros(2*nz, nx*nz), E_j)
    end
    for j = 1:nz
        f[j+nz+3] = λ -> S(λ, nz+j)
        e_j = zeros(nz)
        e_j[j] = 1
        E_j = [spzeros(nz); R(e_j)]
        E_j = E_j * (E_j/nz)'
        A[j+nz+3] =  hvcat((2,2), spzeros(nx*nz,nx*nz), spzeros(nx*nz, 2*nz), spzeros(2*nz, nx*nz), E_j)
    end
    return SPMF_NEP(A,f)
end    


###########################################################
# Generate R-matrix
# Part of defining the P-matrix, see above, Jarlebring-(1.6) and Ringh-(2.8) and Remark 1
function generate_R_matrix(nz::Integer)
    # The scaled FFT-matrix R
    p = (nz-1)/2;
    bb = exp(-2im*pi*((1:nz)-1)*(-p)/nz);  # scaling to do after FFT
    function R(X)
        return flipdim(bb .* fft(X), 1);
    end
    bbinv = 1./bb; # scaling to do before inverse FFT
    function Rinv(X)
        return ifft(bbinv .* flipdim(X,1));
    end
    return R, Rinv
end


###########################################################
# Generate S-function for matrix argument
# Part of defining the P-matrix, see above, Jarlebring-(2.4) and Ringh-(2.8)(2.3)
function generate_S_function(nz::Integer, hx, Km, Kp)
    # Constants from the problem
    p = (nz-1)/2;
    d0 = -3/(2*hx);
    b = 4*pi*1im * (-p:p);
    cM = Km^2 - 4*pi^2 * ((-p:p).^2);
    cP = Kp^2 - 4*pi^2 * ((-p:p).^2);

    betaM = function(γ, j::Integer)
        return γ^2 + b[j]*γ + cM[j]*eye(Complex128,size(γ,1))
    end
    betaP = function(γ, j::Integer)
        return γ^2 + b[j]*γ + cP[j]*eye(Complex128,size(γ,1))
    end

    sM = function(γ, j::Integer)
        return  1im*sqrtm_schur_pos_imag(betaM(γ, j)) + d0*eye(Complex128, size(γ,1))
    end
    sP = function(γ, j::Integer)
        return  1im*sqrtm_schur_pos_imag(betaP(γ, j)) + d0*eye(Complex128, size(γ,1))
    end

    
    S = function(γ, j::Integer)
        if j <= nz
            return sM(γ,j)
        elseif j <= 2*nz
            return sP(γ,(j-nz))
        else
            error("The chosen j = ", j, "but the setup nz = ", nz, ". Hence j>2nz which is illegal.")
        end
    end

    return S
end


###########################################################
    """
    sqrtm_schur_pos_imag(A::AbstractMatrix)
 Computes the matrix square root on the 'correct branch',
 that is, with positivt imaginary part. Similar to Schur method
 in Algorithm 6.3 in Higham matrix functions. 
"""
function sqrtm_schur_pos_imag(A::AbstractMatrix)
    n = size(A,1);
    AA = Array{Complex128,2}(A);
    (T, Q, ) = schur(AA)
    U = zeros(Complex128,n,n);
    for i = 1:n
        U[i,i] = sqrt_pos_imag(T[i,i])
    end
    private_inner_loops_sqrt!(n, U, T)
    return Q*U*Q'
end
#Helper function executing the inner loop (more Juliaesque)
function private_inner_loops_sqrt!(n, U, T)
    temp = zero(Complex128);
    for j = 2:n
        for i = (j-1):-1:1
            temp *= zero(Complex128);
            for k = (i+1):(j-1)
                temp += U[i,k]*U[k,j]
            end
            U[i,j] = (T[i,j] - temp)/(U[i,i]+U[j,j])
        end
    end
end

    """
    sqrt_pos_imag(a::Complex128) and sqrt_pos_imag(a::Float64)
 Helper function: Computes the scalar square root on the 'correct branch',
 that is, with positivt imaginary part.
"""
function sqrt_pos_imag(a::Complex128)
    return sign(imag(a))*sqrt(a)
end
function sqrt_pos_imag(a::Float64)
    return sqrt(a)
end


###########################################################
# Waveguide eigenvalue problem - WEP
# A more optimized (native) implementation of the WEP with FD discretization
    """
    Waveguide eigenvalue problem
  A more optimized implementation of the WEP for FD-discretization.\\
  Closer to what is done in the article:
    ''E. Ringh, and G. Mele, and J. Karlsson, and E. Jarlebring,
      Sylvester-based preconditioning for the waveguide eigenvalue problem,
      Linear Algebra and its Applications''
"""
    type WEP_FD <: NEP
        nx::Integer
        nz::Integer
        hx
        hz
        A::Function
        B::Function
        C1
        C2T
        k_bar
        K
        p::Integer
        d0
        d1
        d2
        b
        cMP # cM followed by cP in a vector (length = 2*nz)
        R::Function
        Rinv::Function
        Pinv::Function
        generate_Pm_and_Pp_inverses::Function

        function WEP_FD(nx, nz, hx, hz, Dxx, Dzz, Dz, C1, C2T, K, Km, Kp)
            n = nx*nz + 2*nz
            k_bar = mean(K)
            K_scaled = K-k_bar*ones(Complex128,nz,nx)

            A = function(λ, d=0)
                if(d == 0)
                    return Dzz + 2*λ*Dz + λ^2*speye(Complex128, nz, nz) + k_bar*speye(Complex128, nz, nz)
                elseif(d == 1)
                    return 2*Dz + 2*λ*speye(Complex128, nz, nz)
                elseif(d == 2)
                    return 2*speye(Complex128, nz, nz)
                else
                    return spzeros(Complex128, nz, nz)
                end
            end
            B = function(λ, d=0)
                if(d == 0)
                    return Dxx
                else
                    return spzeros(Complex128, nx, nx)
                end
            end

            p = (nz-1)/2

            d0 = -3/(2*hx)
            d1 = 2/hx
            d2 = -1/(2*hx)

            b = 4*pi*1im * (-p:p)
            cM = Km^2 - 4*pi^2 * ((-p:p).^2)
            cP = Kp^2 - 4*pi^2 * ((-p:p).^2)
            cMP = [cM; cP]

            R, Rinv = generate_R_matrix(nz)
            Pinv = generate_Pinv_matrix(nz, hx, Km, Kp)
            generate_Pm_and_Pp_inverses(σ) =  helper_generate_Pm_and_Pp_inverses(nz, b, cMP, d0, R, Rinv, σ)

            this = new(nx, nz, hx, hz, A, B, C1, C2T, k_bar, K_scaled, p, d0, d1, d2, b, cMP, R, Rinv, Pinv, generate_Pm_and_Pp_inverses)
        end
    end


    function size(nep::WEP_FD, dim=-1)
        n = nep.nx*nep.nz + 2*nep.nz
        if (dim==-1)
            return (n,n)
        else
            return n
        end
    end


    function issparse(nep::WEP_FD)
        return false
    end


    #Helper function: Generates function to compute iverse of the boundary operators, Ringh - (2.8)
    #To be used in the Schur-complement- and SMW-context.
    function helper_generate_Pm_and_Pp_inverses(nz, b, cMP, d0, R, Rinv, σ)
        # S_k(σ) + d_0, as in Ringh - (2.3a)
        coeffs = zeros(Complex128, 2*nz)
            aa = 1.0
        for j = 1:2*nz
            bb = b[rem(j-1,nz)+1]
            cc = cMP[j]
            coeffs[j] = 1im*sqrt_derivative(aa, bb, cc, 0, σ) + d0
        end

        # P_inv_m and P_inv_p, the boundary operators
        P_inv_m = function(v)
            return R(Rinv(v) ./ coeffs[1:nz])
        end
        P_inv_p = function(v)
            return R(Rinv(v) ./ coeffs[(nz+1):(2*nz)])
        end

        return P_inv_m, P_inv_p
    end


    """
    compute_Mlincomb(nep::WEP_FD, λ::Number, V; a=ones(Complex128,size(V,2)))
Specialized for Waveguide Eigenvalue Problem discretized with Finite Difference\\\\
 Computes the linear combination of derivatives\\
 ``Σ_i a_i M^{(i)}(λ) v_i``
"""
    function compute_Mlincomb(nep::WEP_FD, λ::Number, V; a=ones(Complex128,size(V,2)))
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
        y1_mat = zeros(Complex128, nz, nx)
        y1_mat += nep.A(λ) * V1_mat[:,:,1]  +  V1_mat[:,:,1] * nep.B(λ)  +  nep.K .* V1_mat[:,:,1]
        for d = 1:min(max_d,3)
            y1_mat += nep.A(λ,d) * V1_mat[:,:,d+1]
        end
        y1 = y1_mat[:]
        y1 += nep.C1 * V2[:,1]

        # Compute the bottom part (2*nz)
        y2 = zeros(Complex128, 2*nz, 1) #Make a (2*nz,1) since application of R returns on that format
        D = zeros(Complex128, 2*nz, na)
        for j = 1:2*nz
            a = 1
            b = nep.b[rem(j-1,nz)+1]
            c = nep.cMP[j]
            der_coeff = 1im*sqrt_derivative(a, b, c, max_d, λ)
            for jj = 1:na
                D[j, jj] = der_coeff[jj]
            end
        end

        y2 += (D[:,1] + nep.d0) .* [nep.Rinv(V2[1:nz,1]); nep.Rinv(V2[nz+1:2*nz,1])] #Multpilication with diagonal matrix optimized by working "elementwise" Jarlebring-(4.6)
        for jj = 2:na
            y2 += D[:,jj] .* [nep.Rinv(V2[1:nz,jj]); nep.Rinv(V2[nz+1:2*nz,jj])] #Multpilication with diagonal matrix optimized by working "elementwise" Jarlebring-(4.6)
        end
        y2 = vec([nep.R(y2[1:nz]); nep.R(y2[nz+1:2*nz])]) #Make a vector, since in Julia the types (2*nz,) and (2*nz,1) are different
        y2 += nep.C2T * V1[:,1] #Action of C2T. OBS: Add last because of implcit storage in R*D_i*R^{-1}*v_i

        return [y1;y2]
    end



###########################################################
#Special instances of function calls for Mlincomb_matvec and GMRES to take into account that
#the linear system to solve is with SCHUR COMPLEMENT. Ringh - Algorithm 2, step 10

    function *{T_num}(M::Mlincomb_matvec{T_num, WEP_FD}, v::AbstractVector)
    # Ringh - (2.13)(3.3)
        λ = M.λ
        nep = M.nep

        P_inv_m, P_inv_p = nep.generate_Pm_and_Pp_inverses(λ)

        X = reshape(v, nep.nz, nep.nx)

        return vec(  vec( nep.A(λ)*X + X*nep.B(λ) + nep.K.*X ) - nep.C1 * nep.Pinv(λ, nep.C2T*v)  )

    end

    function size{T_num}(M::Mlincomb_matvec{T_num, WEP_FD}, dim=-1)
        return M.nep.nx * M.nep.nz
    end


    function lin_solve{T_num}(solver::GMRESLinSolver{T_num, WEP_FD}, x::Array; tol=eps(real(T_num)))
    # Ringh - Proposition 2.1
        λ = solver.A.λ
        nep = solver.A.nep

        x_int = x[1:(nep.nx*nep.nz)]
        x_ext = x[((nep.nx*nep.nz)+1):((nep.nx*nep.nz) + 2*nep.nz)]

        rhs = vec(  x_int - nep.C1*nep.Pinv(λ, x_ext))

        if( solver.gmres_log )
            q, convhist = gmres(solver.A, rhs; tol=tol, solver.kwargs...)
        else
            q = gmres(solver.A, rhs; tol=tol, solver.kwargs...)
        end

        x = [q; vec(nep.Pinv(λ, nep.C2T * q + x_ext))]

        return x
    end

# Generate P^{-1}-matrix
# P is the lower right part of the system matrix, from the DtN maps Jarlebring-(1.5)(1.6) and Ringh-(2.4)(2.8)
function generate_Pinv_matrix(nz::Integer, hx, Km, Kp)

    R, Rinv = generate_R_matrix(nz::Integer)
    p = (nz-1)/2;

    # Constants from the problem
    d0 = -3/(2*hx);
    a = ones(Complex128,nz);
    b = 4*pi*1im * (-p:p);
    cM = Km^2 - 4*pi^2 * ((-p:p).^2);
    cP = Kp^2 - 4*pi^2 * ((-p:p).^2);


    function betaM(γ)
        return a*γ^2 + b*γ + cM
    end
    function betaP(γ)
        return a*γ^2 + b*γ + cP
    end

    function sM(γ::Number)
        bbeta = betaM(γ)
        return 1im*sign(imag(bbeta)).*sqrt(bbeta)+d0;
    end
    function sP(γ::Number)
        bbeta = betaP(γ)
        return 1im*sign(imag(bbeta)).*sqrt(bbeta)+d0;
    end

    # BUILD THE INVERSE OF THE FOURTH BLOCK P
    function P(γ,x::Union{Array{Complex128,1}, Array{Float64,1}})
        return vec(  [R(Rinv(x[1:Int64(end/2)]) ./ sM(γ));
                      R(Rinv(x[Int64(end/2)+1:end]) ./ sP(γ))  ]  )
    end


    return P
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

    derivatives = zeros(Complex128,d+1)

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

    yip2 = zero(Complex128)
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

include("waveguide_preconditioner.jl")

end
