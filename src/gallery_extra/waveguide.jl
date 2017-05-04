

export gallery_waveguide
export matlab_debug_WEP_FD #ONLY FOR DEBUGGING
export matlab_debug_full_matrix_WEP_FD_SPMF #ONLY FOR DEBUGGING
export debug_sqrtm_schur #ONLY FOR DEBUGGING
export fft_debug_mateq #ONLY FOR DEBUGGING
export debug_sqrt_derivative #ONLY FOR DEBUGGING

# Specializalized NEPs
export WEP_FD



 # We overload these
    import NEPCore.compute_Mder
    import NEPCore.compute_Mlincomb
#    import NEPCore.compute_MM
#    import NEPCore.compute_resnorm
#    import NEPCore.compute_rf

    import Base.size
    import Base.issparse
    import Base.*
    import Base.norm

    export compute_Mder
    export compute_Mlincomb
#    export compute_MM
#    export compute_resnorm
#    export compute_rf
    export size
    export issparse
    export *
    export norm




"
Creates the NEP associated with example in
E. Ringh, and G. Mele, and J. Karlsson, and E. Jarlebring, 
Sylvester-based preconditioning for the waveguide eigenvalue problem,
Linear Algebra and its Applications, 2017

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
        K, hx, hz, k = generate_wavenumber_fd( nx, nz, waveguide, delta)
        P, p_P = generate_P_matrix(nz, hx, k)
        Dxx, Dzz, Dz = generate_fd_interion_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        nep = assemble_waveguide_spmf_fd(nx, nz, hx, Dxx, Dzz, Dz, C1, C2T, K, k)

    elseif (NEP_format == "WEP") && (discretization == "FD")
        K, hx, hz, k = generate_wavenumber_fd( nx, nz, waveguide, delta)
        P, p_P = generate_P_matrix(nz, hx, k)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        nep = WEP_FD

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
function assemble_waveguide_spmf_fd(nx::Integer, nz::Integer, hx, Dxx::SparseMatrixCSC, Dzz::SparseMatrixCSC, Dz::SparseMatrixCSC, C1::SparseMatrixCSC, C2T::SparseMatrixCSC, K::Union{Array{Complex128,2},Array{Float64,2}}, k::Function)
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
    f[1] = λ -> 1
    f[2] = λ -> λ
    f[3] = λ -> λ^2

    R, Rinv = generate_R_matrix(nz)
    S = generate_S_function(nz, hx, k)

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
# Generate discretization matrices for FINITE DIFFERENCE
    """
 Genearate the discretization matrices for the interior for  Finite Difference.
"""
function generate_fd_interion_mat( nx, nz, hx, hz)
    ex = ones(nx)
    ez = ones(nz)

    # DISCRETIZATION OF THE SECOND DERIVATIVE
    Dxx = spdiagm((ex[1:end-1], -2*ex, ex[1:end-1]), (-1, 0, 1), nx, nx)
    Dzz = spdiagm((ez[1:end-1], -2*ez, ez[1:end-1]), (-1, 0, 1), nz, nz)
    #IMPOSE PERIODICITY IN Z-DIRECTION
    Dzz[1, end] = 1;
    Dzz[end, 1] = 1;

    Dxx = Dxx/(hx^2);
    Dzz = Dzz/(hz^2);

    # DISCRETIZATION OF THE FIRST DERIVATIVE
    Dz  = spdiagm((-ez[1:end-1], ez[1:end-1]), (-1, 1), nz, nz);

    #IMPOSE PERIODICITY
    Dz[1, end] = -1;
    Dz[end, 1] = 1;

    Dz = Dz/(2*hz);

    return Dxx, Dzz, Dz
end
    """
 Genearate the discretization matrices for Finite Difference.
"""
function generate_fd_boundary_mat( nx, nz, hx, hz)
    # BUILD THE SECOND BLOCK C1
    e1 = spzeros(nx,1)
    e1[1] = 1
    en = spzeros(nx,1)
    en[end] = 1
    Iz = speye(nz,nz)
    C1 = [kron(e1,Iz) kron(en,Iz)]/(hx^2);


    # BUILD THE THIRD BLOCK C2^T
    d1 = 2/hx;
    d2 = -1/(2*hx);
    vm = spzeros(1,nx);
    vm[1] = d1;
    vm[2] = d2;
    vp = spzeros(1,nx);
    vp[end] = d1;
    vp[end-1] = d2;
    C2T = [kron(vm,Iz); kron(vp,Iz)];

    return C1, C2T
end


###########################################################
# Generate Wavenumber FINITE DIFFERENCE
    """
 Genearate a wavenumber for Finite Difference.
"""
function generate_wavenumber_fd( nx::Integer, nz::Integer, wg::String, delta::Number)
    if wg == "TAUSCH"
        return generate_wavenumber_fd_tausch( nx, nz, delta)
    elseif wg == "JARLEBRING"
        return generate_wavenumber_fd_jarlebring( nx, nz, delta) 
    end
    # Use early-bailout principle. If a supported waveguide is found, compute and return. Otherwise end up here and throw and error
    error("No wavenumber loaded: The given Waveguide '", wg ,"' is not supported in 'FD' discretization.")
end


    """
 Genearate the wavenumber in FD discretization for the waveguide
described by TAUSCH.
"""
function generate_wavenumber_fd_tausch( nx::Integer, nz::Integer, delta::Number)
    xm = 0;
    xp = (2/pi) + 0.4;
    zm = 0;
    zp = 1;

    xm = xm - delta;
    xp = xp + delta;

    # Domain (First generate including the boundary)
    X = linspace(xm, xp, nx+2);
    const hx = step(X);
    X = collect(X);
    Z =linspace(zm, zp, nz+1);
    const hz = step(Z);
    Z = collect(Z);
    # Removing the boundary
    X = X[2:end-1];
    Z = Z[2:end];


    # The actual wavenumber
    const k1 = sqrt(2.3)*pi;
    const k2 = sqrt(3)*pi;
    const k3 = pi;
    k = function(x,z)
            z_ones = ones(size(z)) #Placeholder of ones to expand x-vector
            k1*(x .<= 0) .* z_ones +
            k2*(x.>0) .* (x.<=2/pi) .* z_ones +
            k2*(x.>2/pi) .* (x.<=(2/pi+0.4)) .* (z.>0.5) +
            k3*(x.>2/pi) .* (z.<=0.5) .* (x.<=(2/pi+0.4)) +
            k3*(x.>(2/pi+0.4)) .* z_ones;
        end

    const K = k(X', Z).^2;

    return K, hx, hz, k
end


    """
 Genearate the wavenumber in FD discretization for the waveguide
described by JARLEBRING.
"""
function generate_wavenumber_fd_jarlebring( nx::Integer, nz::Integer, delta::Number)
    xm = -1;
    xp = 1;
    zm = 0;
    zp = 1;

    xm = xm - delta;
    xp = xp + delta;

    # Domain (First generate including the boundary)
    X = linspace(xm, xp, nx+2);
    const hx = step(X);
    X = collect(X);
    Z =linspace(zm, zp, nz+1);
    const hz = step(Z);
    Z = collect(Z);
    # Removing the boundary
    X = X[2:end-1];
    Z = Z[2:end];


    # The actual wavenumber
    const k1 = sqrt(2.3)*pi;
    const k2 = 2*sqrt(3)*pi;
    const k3 = 4*sqrt(3)*pi;
    const k4 = pi;
    k = function(x,z)
            z_ones = ones(size(z)) #Placeholder of ones to expand x-vector
            x_ones = ones(size(x)) #Placeholder of ones to expand z-vector
            k1 *(x.<=-1) .* z_ones  +
            k4 *(x.>1) .* z_ones  +
            k4 *(x.>(-1+1.5)) .* (x.<=1) .* (z.<=0.4) +
            k3 *(x.>(-1+1)) .* (x.<=(-1+1.5)) .* z_ones +
            k3 *(x.>(-1+1.5)) .* (x.<=1) .* (z.>0.4) +
            k3 *(x.>-1) .* (x.<=(-1+1)) .* (z.>0.5) .* (z.*x_ones-(x.*z_ones)/2.<=1) +
            k2 *(x.>-1) .* (x.<=(-1+1)) .* (z.>0.5) .* (z.*x_ones-(x.*z_ones)/2.>1) +
            k3 *(x.>-1) .* (x.<=(-1+1)) .* (z.<=0.5) .* (z.*x_ones+(x.*z_ones)/2.>0) +
            k2 *(x.>-1) .* (x.<=(-1+1)) .* (z.<=0.5) .* (z.*x_ones+(x.*z_ones)/2.<=0);
        end

    const K = k(X', Z).^2;
    return K, hx, hz, k
end


###########################################################
# Generate discretization matrices for FINITE ELEMENT
    """
 Genearate the discretization matrices for Finite Difference.
"""
function generate_fem_matrices( nx, nz, hx, hz)
        error("FEM discretization currently not supported.")
end


###########################################################
# Generate Wavenumber FINITE ELEMENT
    """
 Genearate a wavenumber for Finite Element.
"""
function generate_wavenumber_fem( nx::Integer, nz::Integer, wg::String, delta::Number)
    # Use early-bailout principle. If a supported waveguide is found, compute and return. Otherwise end up here and throw and error (see FD implementation)
    error("No wavenumber loaded: The given Waveguide '", wg ,"' is not supported in 'FEM' discretization.")
end


###########################################################
# Generate P-matrix
function generate_P_matrix(nz::Integer, hx, k::Function)

    R, Rinv = generate_R_matrix(nz::Integer)
    const p = (nz-1)/2;

    # Constants from the problem
    const Km = k(-Inf, 1/2)[1];
    const Kp = k(Inf, 1/2)[1];
    const d0 = -3/(2*hx);
    const a = ones(Complex128,nz);
    const b = 4*pi*1im * (-p:p);
    const cM = Km^2 - 4*pi^2 * ((-p:p).^2);
    const cP = Kp^2 - 4*pi^2 * ((-p:p).^2);


    function betaM(γ)
        return a*γ^2 + b*γ + cM
    end
    function betaP(γ)
        return a*γ^2 + b*γ + cP
    end

    const signM = 1im*sign(imag(betaM(-1-1im))); # OBS! LEFT HALF-PLANE!
    const signP = 1im*sign(imag(betaP(-1-1im))); # OBS! LEFT HALF-PLANE!

    function sM(γ::Number)
        return signM.*sqrt(betaM(γ))+d0;
    end
    function sP(γ::Number)
        return signP.*sqrt(betaP(γ))+d0;
    end

    function p_sM(γ)
        return signM.*(2*a*γ+b)./(2*sqrt(a*γ^2+b*γ+cM));
    end
    function p_sP(γ)
        return signP.*(2*a*γ+b)./(2*sqrt(a*γ^2+b*γ+cP));
    end

    # BUILD THE FOURTH BLOCK P
    function P(γ,x::Union{Array{Complex128,1}, Array{Float64,1}})
        return [R(Rinv(x[1:Int64(end/2)]) .* sM(γ));
                R(Rinv(x[Int64(end/2)+1:end]) .* sP(γ))  ];
    end

    # BUILD THE DERIVATIVE OF P
    function p_P(γ,x::Union{Array{Complex128,1}, Array{Float64,1}})
        return [R(Rinv(x[1:Int64(end/2)]) .* p_sM(γ));
                R(Rinv(x[Int64(end/2)+1:end]) .* p_sP(γ))  ];
    end

    return P, p_P
end


###########################################################
# Generate R-matrix
function generate_R_matrix(nz::Integer)
    # The scaled FFT-matrix R
    const p = (nz-1)/2;
    const bb = exp(-2im*pi*((1:nz)-1)*(-p)/nz);  # scaling to do after FFT
    function R(X)
        return flipdim((bb*ones(size(X,2),1)') .* fft(X), 1);
    end
    bbinv = 1./bb; # scaling to do before inverse FFT
    function Rinv(X)
        return ifft((bbinv*ones(size(X,2),1)') .* flipdim(X,1));
    end
    return R, Rinv
end


###########################################################
# Generate S-function for matrix argument
function generate_S_function(nz::Integer, hx, k::Function)
    # Constants from the problem
    const p = (nz-1)/2;
    const Km = k(-Inf, 1/2)[1];
    const Kp = k(Inf, 1/2)[1];
    const d0 = -3/(2*hx);
    const b = 4*pi*1im * (-p:p);
    const cM = Km^2 - 4*pi^2 * ((-p:p).^2);
    const cP = Kp^2 - 4*pi^2 * ((-p:p).^2);

    function betaM(γ::AbstractArray, j::Integer)
        return γ^2 + b[j]*γ + cM[j]*speye(size(γ,1))
    end
    function betaP(γ::AbstractArray, j::Integer)
        γ^2 + b[j]*γ + cP[j]*speye(size(γ,1))
    end

    function sM(γ::AbstractArray, j::Integer)
        return  1im*speye(Complex128, size(γ,1)) * sqrtm_schur_pos_imag(betaM(γ, j)) + d0*speye(Complex128, size(γ,1))
    end
    function sP(γ::AbstractArray, j::Integer)
        return  1im*speye(Complex128, size(γ,1)) * sqrtm_schur_pos_imag(betaP(γ, j)) + d0*speye(Complex128, size(γ,1))
    end

    function S(γ::AbstractArray, j::Integer)
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
# Compute the matrix square root on the "correct branch",
# that is, with positive imaginary part
function sqrtm_schur_pos_imag(A::AbstractMatrix)
    n = size(A,1);
    AA = full(complex(A))
    (T, Q, ) = schur(AA)
    U = zeros(Complex128,n,n);
    for i = 1:n
        U[i,i] = sign(imag(T[i,i]))*sqrt(T[i,i])
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




###########################################################
# Waveguide eigenvalue problem - WEP
# A more optimized (native) implementation of the WEP with FD discretization

    """
### Waveguide eigenvalue problem
  A more optimized implementation of the WEP for FD-discretization.\\
  Closer to what is done in the article:
    ''E. Ringh, and G. Mele, and J. Karlsson, and E. Jarlebring,
      Sylvester-based preconditioning for the waveguide eigenvalue problem,
      Linear Algebra and its Applications, 2017''
"""
    type WEP_FD <: NEP
        nx::Integer
        nz::Integer
    end


    function size(nep::WEP_FD, dim=-1)
        n = nep.nx * nep.nz + 2*nep.nz
        if (dim==-1)
            return (n,n)
        else
            return n
        end
    end


    function issparse(nep::WEP_FD)
        return false
    end


    """
    compute_Mder(nep::WEP_FD, λ::Number, i::Integer=0)
 Gives an object that acts as a matrix-multiplication
 for the WEP.
 Higher order derivatives gives an error.
"""
     # TODO: This function compute only the 0:th and 1:st derivatives. Extend?
    function compute_Mder(nep::WEP_FD, λ::Number, i::Integer=0)
        if(i == 0)
            return compute_WEP_matrix_object(nep, λ)
        elseif( i == 1)
            
        else
            error("Cannot compute derivative for higher order of WEP.")
        end
    end


    function compute_WEP_matrix_object(nep::WEP_FD, λ::Number)
        
    end

    """
    An abstract matrix object from the WEP_FD.\\
    Overload * and size() to make it act like a normal matrix
"""
    type WEP_matrix_object
        
    end

    function *(M::WEP_matrix_object, v::AbstractVector)
        
    end

    function size(M::WEP_matrix_object)
        
    end

    function norm(M::WEP_matrix_object, p::Real=1)
        
    end


    """
    compute_Mlincomb(nep::WEP_FD, λ::Number, V; a=ones(Complex128,size(V,2)))
Specialized for Waveguide Eigenvalue Problem discretized with Finite Difference\\\\
 Computes the linear combination of derivatives\\
 ``Σ_i a_i M^{(i)}(λ) v_i``
"""
    function compute_Mlincomb(nep::WEP_FD, λ::Number, V; a=ones(Complex128,size(V,2)))
        na = size(a)
        nv, mv = size(V)
        n_nep = size(nep,1)
        if(na != mv)
            error("Incompatible sizes: Number of coefficients = ", na, ", number of vectors = ", mv, ".")
        end
        if(nv != n_nep)
            error("Incompatible sizes: Length of vectors = ", nv, ", size of NEP = ", n_nep, ".")
        end
    end



    """
    fft_wg( C, gamma, kk, hx, hz )
 Solves the Sylvester equation for the WEP.
"""
function solve_wg_sylvester_fft( C, λ, kk, hx, hz )

    nz = size(C,1)
    nx = size(C,2)

    alpha = λ^2+kk;

    #eigenvalues of A
    v=zeros(nz);   v[1]=-2;    v[2]=1;     v[nz]=1;   v=v/(hz^2);
    w=zeros(nz);   w[2]=1;     w[nz]=-1;   w=w*(λ/hz);
    D=fft(v+w)+alpha;

    # eigenvalues of B = Dxx
    S = -(4/hx^2) * sin(pi*(1:nx)/(2*(nx+1))).^2
#    S=S.'

    # solve the diagonal matrix equation
    Z=zeros(Complex128,nz,nx)

    CC = Vh!( Wh(C')' )

    for k=1:nx
        Z[:,k] += CC[:,k]./(D+S[k])
    end


    # change variables
    return V!((W(Z'))')

end


# Start: Auxiliary computations of eigenvector actions using FFT
    function V!(X)
    # Compute the action of the eigenvectors of A = Dzz + Dz + c*I
        nx = size(X,2)
        return fft!(X,1)/sqrt(nx)
    end

    function Vh!(X)
    # Compute the action of the transpose of the eigenvectors of A = Dzz + Dz + c*I
        nx = size(X,2)
        return ifft!(X,1)*sqrt(nx)
    end

    function W( X )
    #W Compute the action of the matrix W
    #   W is the matrix of the eigenvectors of the second derivative Dxx
    #   W*X can be computed with FFTs

        WX = (1im/2)*(F(X)-Fh(X))

        nz = size(X,1)
        return WX/sqrt((nz+1)/2)

    end

    function Wh( X )
    #Wh Compute the action of the matrix Wh
    #   Wh is the transpose of the matrix of the eigenvectors of the second derivative Dxx
    #   Wh*X can be computed with FFTs

        WX = (Fh(X)-F(X))/(2im)

        nz = size(X,1)
        return WX/sqrt((nz+1)/2)

    end

    function F( v )
    #F is an auxiliary function for W and Wh

        m=size(v,2)

        v=[zeros(1,m); v]
        n=size(v,1)

        pad = [v; zeros(eltype(v),n,m)]
        v = fft!(pad,1)
        return v[2:n,:]

    end

    function Fh( v )
    #Fh is an auxiliary function for W and Wh

        m=size(v,2)

        v=[zeros(1,m); v]
        n=size(v,1)

        pad = [v; zeros(eltype(v),n,m)]
        v = ifft!(pad,1)
        return v[2:n,:]*2*n

    end
# End: Auxiliary computations of eigenvector actions using FFT






######################## DEBUG ############################
###########################################################
# DEBUG: Test the generated matrices against MATLAB code
using MATLAB
function matlab_debug_WEP_FD(nx::Integer, nz::Integer, delta::Number)
    if(nx > 200 || nz > 200)
        warn("This debug is 'naive' and might be slow for the discretization used.")
    end

    #The error observed likely comes from difference in linspace-implementation.
    #include("../bugs/test_linspace.jl")

    γ = -rand() - 1im*rand()
    gamma = γ


    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing waveguide: ", waveguide)

        K, hx, hz, k = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interion_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        P, p_P = generate_P_matrix(nz, hx, k)

        R, Rinv = generate_R_matrix(nz)
        S = generate_S_function(nz, hx, k)
        P_j2 = zeros(Complex128, 2*nz,2*nz)
        D1 = zeros(Complex128, nz,nz)
        D2 = zeros(Complex128, nz,nz)
        for j = 1:nz
            D1[j,j] = S([γ]'',j)[1]
        end
        for j = 1:nz
            D2[j,j] = S([γ]'',j+nz)[1]
        end
        Iz = eye(nz,nz);
        P_j2[1:nz,1:nz] = R(D1*Rinv(Iz))
        P_j2[(nz+1):(2*nz),(nz+1):(2*nz)] = R(D2*Rinv(Iz))


        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = "../matlab/WEP"
        @mput nx nz delta WEP_path waveguide_str gamma
        @matlab begin
            addpath(WEP_path)
            nxx = double(nx)
            nzz = double(nz)
            options = struct
            options.delta = delta
            options.wg = waveguide_str
            matlab_nep = nep_wg_generator(nxx, nzz, options)

            P_m = NaN(2*nzz, 2*nzz);
            Iz = eye(2*nzz);
            eval("for i = 1:2*nzz;   P_m(:,i) = matlab_nep.P(gamma, Iz(:,i));    end")
            C1_m = matlab_nep.C1;
            C2T_m = matlab_nep.C2T;
            K_m = matlab_nep.K;
            hx_m = matlab_nep.hx;
            hz_m = matlab_nep.hz;

        @matlab end
        @mget K_m C2T_m C1_m hx_m hz_m P_m
        println("  -- Matlab printouts end --")

        P_j = zeros(Complex128, 2*nz,2*nz)
        Iz = eye(2*nz, 2*nz)
        for i = 1:2*nz
            P_j[:,i] = P(γ, Iz[:,i])
        end

        println("Difference hx_m - hx = ", abs(hx_m-hx))
        println("Relative difference (hx_m - hx)/hx = ", abs(hx_m-hx)/abs(hx))
        println("Difference hz_m - hz = ", abs(hz_m-hz))
        println("Difference K_m  -K = ", norm(K_m-K))
        println("Difference C1_m - C1 = ", norm(full(C1_m-C1)))
        println("Relative difference norm(C1_m - C1)/norm(C1) = ", norm(full(C1_m-C1))/norm(full(C1)))
        println("Difference C2T_m - C2T = ", norm(full(C2T_m-C2T)))
        println("Relative difference norm(C2T-m - C2T)/norm(C2T) = ", norm(full(C2T_m-C2T))/norm(full(C2T)))
        println("Difference P_m(γ) - P(γ) = ", norm(P_m-P_j))
        println("Relative difference norm(P_m(γ) - P(γ))/norm(P(γ)) = ", norm(P_m-P_j)/norm(P_j))
        println("Difference P_m(γ) - P_2(γ) = ", norm(P_m-P_j2))
        println("Relative difference norm(P_m(γ) - P_2(γ))/norm(P_2(γ)) = ", norm(P_m-P_j2)/norm(full(P_j2)))
    end
end

# Test the full generated system-matrix against against MATLAB code
function matlab_debug_full_matrix_WEP_FD_SPMF(nx::Integer, nz::Integer, delta::Number)
    if(nx > 40 || nz > 40)
        warn("This debug is 'naive' and might be slow for the discretization used.")
    end


    γ = -rand() - 1im*rand()
    gamma = γ

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing full matrix M for waveguide: ", waveguide)

        nep_j = nep_gallery("waveguide", nx, nz, waveguide, "fD", "SpmF", delta)
        M_j = compute_Mder(nep_j,γ)

        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = "../matlab/WEP"
        @mput nx nz delta WEP_path waveguide_str gamma
        @matlab begin
            addpath(WEP_path)
            nxx = double(nx)
            nzz = double(nz)
            options = struct
            options.delta = delta
            options.wg = waveguide_str
            matlab_nep = nep_wg_generator(nxx, nzz, options)

            Ixz = eye(nxx*nzz+2*nzz);
            M_m = NaN(nxx*nzz+2*nzz, nxx*nzz+2*nzz);
            eval("for i = 1:(nxx*nzz+2*nzz);   M_m(:,i) = matlab_nep.M(gamma, Ixz(:,i));    end")


        @matlab end
        @mget M_m
        println("  -- Matlab printouts end --")

        println("Difference M_m(γ) - M(γ) = ", norm(full(M_m-M_j)))
        println("Relative difference norm(M_m(γ) - M(γ))/norm(M(γ)) = ", norm(full(M_m-M_j))/norm(full(M_j)))
    end
end


# Test the full generated system-matrix against against MATLAB code
function fft_debug_mateq(nx::Integer, nz::Integer, delta::Number)

    γ = -rand(Complex128)
    gamma = γ
    C = rand(Complex128, nz, nx);
    waveguide = "JARLEBRING"


    K, hx, hz, k = generate_wavenumber_fd( nx, nz, waveguide, delta)
    Dxx, Dzz, Dz = generate_fd_interion_mat( nx, nz, hx, hz)

    k_bar = mean(K)

    A = full(Dzz + 2*γ*Dz + (γ^2+k_bar)*speye(Complex128, nz,nz))
    B = complex(full(Dxx))

    println("Built-in Sylvester solver")
    XX = @time sylvester(A,B,-C)
    println("Relative residual norm = ", norm(A*XX+XX*B-C)/norm(C))

    println("FFT-based Sylvester solver for WG")
    X = @time solve_wg_sylvester_fft( C, γ, k_bar, hx, hz )
    println("Relative residual norm = ", norm(A*X+X*B-C)/norm(C))

end


###########################################################
# Compute the matrix square root
# (only reference implementation, see sqrtm_schur_pos_imag)
function sqrtm_schur(A::AbstractMatrix)
    n = size(A,1);
    (T, Q) = schur(complex(A))
    U = zeros(Complex128,n,n);
    for i = 1:n
        U[i,i] = sqrt(T[i,i])
    end
    for j = 2:n
        for i = (j-1):-1:1
            temp = zero(Complex128)
            for k = (i+1):(j-1)
                temp += U[i,k]*U[k,j]
            end
            U[i,j] = (T[i,j] - temp)/(U[i,i]+U[j,j])
        end
    end
    return Q*U*Q'
end


function debug_sqrtm_schur(n::Integer)
    println("\n\n--- Debugging square root implementations ---\n")
    A = rand(n,n) + 0.1im*rand(n,n);
    sqrtA = sqrtm(A);
    sqrtA2 = sqrtm_schur(A);

    println("Relative error between sqrtm and Schur-fact: ", norm(sqrtA-sqrtA2)/norm(sqrtA))
    println("Relative error between Schur-fact² and A: ", norm(A - sqrtA2^2)/norm(A))

    sqrtA3 = sqrtm_schur_pos_imag(A);
    println("Relative error between Schur-fact-pos-imag² and A: ", norm(A - sqrtA3^2)/norm(A))
    (v,) = eig(sqrtA3)
    test_var = zeros(n)
    TOL = 1e-15;
    for i = 1:n
        test_var[i] = (sign(imag(v[i])) > 0) || abs(imag(v[i])) < TOL
    end
    if sum(test_var) == n
        println("All eigenvalues have negative real part or absolut value smaller than ", TOL)
    else
        println("Eigenvalues with negative real part and absolut value larger than ", TOL, " :")
        for i = 1:n
            if (test_var[i] == 0)
                println("   ", v[i])
            end
        end
    end
    println("\n--- End square root implementations ---\n")
end

###########################################################
# Compute derivative <d> of sqrt(ax^2 + bx + c) in <x>
# (only reference implementation, compute_mlincomb)
function sqrt_derivative(a,b,c, d, x)
    if(d<0)
        error("Cannot take negative derivative. d = ", d)
    end

    aa = a
    bb = b + 2*a*x
    cc = c + a*x^2 + b*x

    yi = sqrt(cc)
    if( d==0 )
        return yi
    end

    yip1 = bb/(2*sqrt(cc))
    fact = Float64(1)
    if( d==1 )
        return yip1 * fact
    end

    yip2 = zero(Complex128)
    for i = 2:d
        m = i - 2
        yip2 = - (2*aa*(m-1)*yi  +  bb*(1+2*m)*yip1) / (2*cc*(2+m))
        fact *= i

        yi = yip1
        yip1 = yip2
    end
    return yip2 * fact
end

function debug_sqrt_derivative()
    println("\n\n--- Debugging derivatives of square root of polynomials ---\n")
    a = rand()
    b = rand()
    c = rand()
    d_vec = [0 1 2 3 4 11 19 20 21 22 30 35 45 60] #Factorial for Int64 overflows at 21!
    x = 25*rand()

        WEP_path = "../matlab/WEP"
        println("  -- Matlab printouts start --")
        @mput a b c d_vec x WEP_path
        @matlab begin
            addpath(WEP_path)
            der_val = sqrt_derivative_test(a,b,c, d_vec, x);
        @matlab end
        @mget der_val
        println("  -- Matlab printouts end --")

    for i = 1:size(d_vec,2)
        d = d_vec[i]
        println("Derivative number d = ", d)
        println("  MATLAB symbolic = ", der_val[i])
        julia_der = sqrt_derivative(a,b,c, d, x)
        println("  Implemented recursion = ", julia_der)
        println("  Relative error = ", abs(der_val[i]-julia_der)/abs(julia_der))
    end

    println("\n--- End derivatives of square root of polynomials ---\n")
end
