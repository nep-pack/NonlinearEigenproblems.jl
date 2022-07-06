"""
    waveguide_debug

A module with an implementation of the waveguide functionality
using MATLAB callbacks.
Used as rerefence implementation to check that Julia implementation works.
"""
module waveguide_debug
push!(LOAD_PATH, string(@__DIR__,"/../"))	# looks for modules in the current directory
push!(LOAD_PATH, string(@__DIR__,"/../../"))	# looks for modules in the current directory

using MATLAB

using Statistics
using Random
using LinearAlgebra
using NonlinearEigenproblems
using NonlinearEigenproblems.NEPCore
using NonlinearEigenproblems.NEPTypes
using NonlinearEigenproblems.LinSolvers
using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using GalleryWaveguide


export matlab_debug_WEP_FD
export matlab_debug_full_matrix_WEP_FD_SPMF
export debug_sqrt_schur
export fft_debug_mateq
export debug_sqrt_derivative
export matlab_debug_Schur_WEP_FD
export debug_Mlincomb_FD_WEP
export debug_Sylvester_SMW_WEP
export matlab_debug_eigval_comp_WEP_FD_and_SPMF
export debug_WEP_FD_preconditioner
export debug_eigval_comp_WEP_FD

#OBS: Explicit imports of non-exported functions
import GalleryWaveguide.generate_wavenumber_fd
import GalleryWaveguide.generate_fd_interior_mat
import GalleryWaveguide.generate_fd_boundary_mat
import GalleryWaveguide.generate_R_matvecs
import GalleryWaveguide.generate_S_function
import GalleryWaveguide.sqrt_schur_pos_imag
import GalleryWaveguide.sqrt_derivative
import GalleryWaveguide.SchurMatVec
import GalleryWaveguide.solve_wg_sylvester_fft!
import GalleryWaveguide.generate_smw_matrix
import GalleryWaveguide.solve_smw
import GalleryWaveguide.construct_WEP_schur_complement


########### SOME REFERENCE IMPLEMENTATIONS ################
###########################################################

# Generate P-matrix
# Is the lower right part of the system matrix, from the DtN maps Jarlebring-(1.5)(1.6) and Ringh-(2.4)(2.8)
function generate_P_matrix(nz::Integer, hx, Km, Kp)

    R, Rinv = generate_R_matvecs(nz::Integer)
    p = (nz-1)/2;

    # Constants from the problem
    d0 = -3/(2*hx);
    a = ones(ComplexF64,nz);
    b = 4*pi*1im * (-p:p);
    cM = Km^2 .- 4*pi^2 * ((-p:p).^2);
    cP = Kp^2 .- 4*pi^2 * ((-p:p).^2);


    function betaM(γ)
        return a*γ^2 + b*γ + cM
    end
    function betaP(γ)
        return a*γ^2 + b*γ + cP
    end

    signM = 1im*sign.(imag(betaM(-1-1im))); # OBS! LEFT HALF-PLANE!
    signP = 1im*sign.(imag(betaP(-1-1im))); # OBS! LEFT HALF-PLANE!

    function sM(γ::Number)
        return signM.*sqrt.(betaM(γ)) .+d0;
    end
    function sP(γ::Number)
        return signP.*sqrt.(betaP(γ)) .+d0;
    end

    function p_sM(γ)
        return signM.*(2*a*γ+b)./(2*sqrt(a*γ^2+b*γ+cM));
    end
    function p_sP(γ)
        return signP.*(2*a*γ+b)./(2*sqrt(a*γ^2+b*γ+cP));
    end

    # BUILD THE FOURTH BLOCK P
    function P(γ,x::Union{Vector{ComplexF64}, Vector{Float64}})
        return vec( [R(Rinv(x[1:Int(end/2)]) .* sM(γ));
                     R(Rinv(x[Int(end/2)+1:end]) .* sP(γ))  ])
    end

    # BUILD THE DERIVATIVE OF P
    function p_P(γ,x::Union{Vector{ComplexF64}, Vector{Float64}})
        return vec( [R(Rinv(x[1:Int(end/2)]) .* p_sM(γ));
                     R(Rinv(x[Int(end/2)+1:end]) .* p_sP(γ))  ])
    end

    return P, p_P
end

# Generate a function for mat-vecs with P^{-1}-matrix
# P is the lower right part of the system matrix, from the DtN maps Jarlebring-(1.5)(1.6) and Ringh-(2.4)(2.8)
function generate_Pinv_matrix(nz::Integer, hx, Km, Kp)

   R, Rinv = generate_R_matvecs(nz::Integer)
   p = (nz-1)/2;

   # Constants from the problem
   d0 = -3/(2*hx);
   a = ones(ComplexF64,nz);
   b = 4*pi*1im * (-p:p);
   cM = Km^2 .- 4*pi^2 * ((-p:p).^2);
   cP = Kp^2 .- 4*pi^2 * ((-p:p).^2);

   function betaM(γ)
       return a*γ^2 + b*γ + cM
   end
   function betaP(γ)
       return a*γ^2 + b*γ + cP
   end

   function sM(γ::Number)
       bbeta = betaM(γ)
       return 1im*sign.(imag(bbeta)).*sqrt.(bbeta) .+ d0
   end
   function sP(γ::Number)
       bbeta = betaP(γ)
       return 1im*sign.(imag(bbeta)).*sqrt.(bbeta) .+ d0
   end

   # BUILD THE INVERSE OF THE FOURTH BLOCK P
   function Pinv(γ,x::Union{Vector{ComplexF64}, Vector{Float64}})
       return vec(  [R(Rinv(x[1:Int(end/2)]) ./ sM(γ));
                     R(Rinv(x[Int(end/2)+1:end]) ./ sP(γ))  ]  )
   end


   return Pinv
end


###########################################################
# Compute the matrix square root
# (only reference implementation, see sqrt_schur_pos_imag)
function sqrt_schur(A::AbstractMatrix)
    n = size(A,1);
    (T, Q) = schur(complex(A))
    U = zeros(ComplexF64,n,n);
    for i = 1:n
        U[i,i] = sqrt(T[i,i])
    end
    for j = 2:n
        for i = (j-1):-1:1
            temp = zero(ComplexF64)
            for k = (i+1):(j-1)
                temp += U[i,k]*U[k,j]
            end
            U[i,j] = (T[i,j] - temp)/(U[i,i]+U[j,j])
        end
    end
    return Q*U*Q'
end


######################## DEBUG ############################
###########################################################

# Test the basic generated FD matrices against MATLAB code
function matlab_debug_WEP_FD(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging Matrices FD against MATLAB ---\n")
    if(nx > 200 || nz > 200)
        @warn "This debug is 'naive' and might be slow for the discretization used."
    end

    #The error observed likely comes from difference in linspace-implementation.
    γ = -rand() - 1im*rand()
    gamma = γ


    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing waveguide: ", waveguide)

        K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        P, p_P = generate_P_matrix(nz, hx, Km, Kp)
        Pinv = generate_Pinv_matrix(nz, hx, Km, Kp)
        P_j = zeros(ComplexF64, 2*nz,2*nz)
        Iz = Matrix(1.0*I, 2*nz, 2*nz)
        for i = 1:2*nz
            P_j[:,i] = P(γ, Iz[:,i])
        end

        R, Rinv = generate_R_matvecs(nz)
        S = generate_S_function(nz, hx, Km, Kp)
        P_j2 = zeros(ComplexF64, 2*nz,2*nz)
        D1 = zeros(ComplexF64, nz,nz)
        D2 = zeros(ComplexF64, nz,nz)
        for j = 1:nz
            D1[j,j] = S(reshape([γ],1,1),j)[1]
        end
        for j = 1:nz
            D2[j,j] = S(reshape([γ],1,1),j+nz)[1]
        end
        Iz = Matrix(1.0*I, nz, nz);
        for j = 1:nz
            P_j2[1:nz,j] = R(D1*Rinv(Iz[:,j]))
            P_j2[(nz+1):(2*nz),j+nz] = R(D2*Rinv(Iz[j,:]))
        end


        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = string(@__DIR__, "/WEP-ref")
        mat"""
            addpath($WEP_path)
            nxx = double($nx);
            nzz = double($nz);
            gamma_m = double($gamma);
            options = [];
            options.delta = $delta;
            options.wg = $waveguide_str;
            matlab_nep = nep_wg_generator(nxx, nzz, options);

            $P_m = NaN(2*nzz, 2*nzz);
            Iz = eye(2*nzz);
            for i = 1:2*nzz
              $P_m(:,i) = matlab_nep.P(gamma_m, Iz(:,i));
            end
            $C1_m = matlab_nep.C1;
            $C2T_m = matlab_nep.C2T;
            $K_m = matlab_nep.K;
            $hx_m = matlab_nep.hx;
            $hz_m = matlab_nep.hz;
        """
        println("  -- Matlab printouts end --")

        println("Difference hx_m - hx = ", abs(hx_m-hx))
        println("Relative difference (hx_m - hx)/hx = ", abs(hx_m-hx)/abs(hx))
        println("Difference hz_m - hz = ", abs(hz_m-hz))
        println("Difference K_m  -K = ", opnorm(K_m-K))
        println("Difference C1_m - C1 = ", opnorm(Matrix(C1_m-C1)))
        println("Relative difference opnorm(C1_m - C1)/opnorm(C1) = ", opnorm(Matrix(C1_m-C1))/opnorm(Matrix(C1)))
        println("Difference C2T_m - C2T = ", opnorm(Matrix(C2T_m-C2T)))
        println("Relative difference opnorm(C2T-m - C2T)/opnorm(C2T) = ", opnorm(Matrix(C2T_m-C2T))/opnorm(Matrix(C2T)))
        println("Difference P_m(γ) - P(γ) = ", opnorm(P_m-P_j))
        println("Relative difference opnorm(P_m(γ) - P(γ))/opnorm(P(γ)) = ", opnorm(P_m-P_j)/opnorm(P_j))
        println("Difference P_m(γ) - P_2(γ) = ", opnorm(P_m-P_j2))
        println("Relative difference opnorm(P_m(γ) - P_2(γ))/opnorm(P_2(γ)) = ", opnorm(P_m-P_j2)/opnorm(Matrix(P_j2)))
        v = rand(ComplexF64,2*nz)
        println("Relative difference norm(v - Pinv_j((γ))*P_j(γ)*v)/norm(v) = ", norm(v - Pinv(γ,P(γ,v)) )/norm(v))
    end
    println("\n--- End Matrices FD against MATLAB ---\n")
end


###########################################################
# Test the full generated system-matrix against MATLAB code
function matlab_debug_full_matrix_WEP_FD_SPMF(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging Full Matrix FD SPMF against MATLAB ---\n")
    if(nx > 40 || nz > 40)
        @warn "This debug is 'naive' and might be slow for the discretization used."
    end


    γ = -rand() - 1im*rand()
    gamma = γ

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing full matrix M for waveguide: ", waveguide)

        nep_j = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "SpmF", delta = delta)
        M_j = compute_Mder(nep_j,γ)

        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = string(@__DIR__, "/WEP-ref")
        mat"""
            addpath($WEP_path)
            nxx = double($nx);
            nzz = double($nz);
            options = [];
            options.delta = double($delta);
            options.wg = $waveguide_str;
            matlab_nep = nep_wg_generator(nxx, nzz, options);

            Ixz = eye(nxx*nzz+2*nzz);
            $M_m = NaN(nxx*nzz+2*nzz, nxx*nzz+2*nzz);
            for i = 1:(nxx*nzz+2*nzz)
              $M_m(:,i) = matlab_nep.M($gamma, Ixz(:,i));
            end
        """
        println("  -- Matlab printouts end --")

        println("Difference M_m(γ) - M(γ) = ", opnorm(Matrix(M_m-M_j)))
        println("Relative difference opnorm(M_m(γ) - M(γ))/opnorm(M(γ)) = ", opnorm(Matrix(M_m-M_j))/opnorm(Matrix(M_j)))
    end
    println("\n--- End Full Matrix FD SPMF against MATLAB ---\n")
end


###########################################################
# Test the square-root implementations. Both reference and negative-imaginary-part.
function debug_sqrt_schur(n::Integer)
    println("\n\n--- Debugging square root implementations ---\n")
    A = rand(n,n) + 0.1im*rand(n,n);
    sqrtA = sqrt(A);
    sqrtA2 = sqrt_schur(A);

    println("Relative error between sqrtm and Schur-fact: ", opnorm(sqrtA-sqrtA2)/opnorm(sqrtA))
    println("Relative error between Schur-fact² and A: ", opnorm(A - sqrtA2^2)/opnorm(A))

    sqrtA3 = sqrt_schur_pos_imag(A);
    println("Relative error between Schur-fact-pos-imag² and A: ", opnorm(A - sqrtA3^2)/opnorm(A))
    v,_ = eigen(sqrtA3)
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
# Test computation of derivative <d> of sqrt(ax^2 + bx + c) in <x>
function debug_sqrt_derivative()
    println("\n\n--- Debugging derivatives of square root of polynomials against MATLAB symbolic ---\n")
    a = 2*rand()
    b = 2*pi*rand(ComplexF64)
    c = 1.67*rand()
    d_vec = [0 1 2 3 4 11 19 20 21 22 30 35 45 60] #Factorial for Int64 overflows at 21!
    x = 25*rand()

        WEP_path = string(@__DIR__, "/WEP-ref")
        println("  -- Matlab printouts start --")
        mat"""
            addpath($WEP_path)
            $der_val = sqrt_derivative_test($a,$b,$c, $d_vec, $x);
        """
        println("  -- Matlab printouts end --")

    julia_der_vec = sqrt_derivative(a,b,c, maximum(d_vec), x)
    for i = 1:size(d_vec,2)
        d = d_vec[i]
        println("Derivative number d = ", d)
        println("    MATLAB symbolic = ", der_val[i])
        julia_der = julia_der_vec[d+1]
        println("    Implemented recursion = ", julia_der)
        println("    Relative error = ", abs(der_val[i]-julia_der)/abs(julia_der))
    end

    println("\n--- End derivatives of square root of polynomials against MATLAB symbolic ---\n")
end


###########################################################
# Test the full generated native-WEP system-matrix against the SPMF-matrix (Mlincomb)
function debug_Mlincomb_FD_WEP(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging Mlincomb (Full Matrix FD WEP-native-format against SPMF) ---\n")
    if(nx > 40 || nz > 40)
        @warn "This debug is 'naive' and might be slow for the discretization used."
    end


    γ = -rand() - 1im*rand()
    gamma = γ
    n = nx*nz+2*nz
    V = rand(ComplexF64, n, 4)

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing full matrix M for waveguide: ", waveguide)

        nep_j = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "WeP", delta = delta)
        nep_SPMF = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "SpmF", delta = delta)

        for d = [0 1 2 5]
            M_j = zeros(ComplexF64, n, n)
            for i = 1:n
                V = zeros(ComplexF64, n)
                V[i] = 1
                M_j[:,i] += compute_Mlincomb(nep_j, γ, V, [1], d)
            end

            M_SPMF = zeros(ComplexF64, n, n)
            for i = 1:n
                V = zeros(ComplexF64, n)
                V[i] = 1
                M_SPMF[:,i] += compute_Mlincomb(nep_SPMF, γ, V, [1], d)
            end
            println("Derivative d = ", d)
            println("    Difference M_SPMF(γ) - M(γ) = ", opnorm(Matrix(M_SPMF-M_j)))
            println("    Relative difference opnorm(M_m(γ) - M(γ))/opnorm(M(γ)) = ", opnorm(Matrix(M_SPMF-M_j))/opnorm(Matrix(M_j)))

            v_j = compute_Mlincomb(nep_j, γ, V)
            v_SPMF = compute_Mlincomb(nep_SPMF, γ, V)
            println("    Relative difference on random set of 4 vectors = ", norm(v_j - v_SPMF)/norm(v_j))
        end


    end
    println("\n--- End Mlincomb (Full Matrix FD WEP-native-format against SPMF) ---\n")
end


###########################################################
# Test the Schur-complement of native-WEP against MATLAB code
function matlab_debug_Schur_WEP_FD(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging Schur-complement of native-WEP against MATLAB ---\n")
    if(nx > 50 || nz > 50)
        @warn "This debug is 'naive' and might be slow for the discretization used."
    end


    γ = -rand() - 1im*rand()
    gamma = γ
    n = nx*nz;

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing Schur-complement for waveguide: ", waveguide)

        nep_j = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "WEp", delta = delta)
        Schur_fun = SchurMatVec(nep_j, γ)
        Schur_j = zeros(ComplexF64, n, n)
        for i = 1:n
            V = zeros(ComplexF64, n)
            V[i] = 1
            Schur_j[:,i] += Schur_fun* V
        end

        Schur_jj = Matrix(construct_WEP_schur_complement(nep_j, γ))

        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = string(@__DIR__, "/WEP-ref")
        mat"""
            addpath($WEP_path)
            nxx = double($nx);
            nzz = double($nz);
            options = [];
            options.delta = double($delta);
            options.wg = $waveguide_str;
            matlab_nep = nep_wg_generator(nxx, nzz, options);

            Ixz = eye(nxx*nzz);
            $Schur_m = NaN(nxx*nzz, nxx*nzz);
            for i = 1:(nxx*nzz)
              $Schur_m(:,i) = matlab_nep.S($gamma, Ixz(:,i));
            end
         """
        println("  -- Matlab printouts end --")

        norm_diff = opnorm(Schur_m-Schur_j)
        println("Difference Schur_m(γ) - Schur_j(γ) = ", norm_diff)
        println("Relative difference opnorm(Schur_m(γ) - Schur_j(γ))/opnorm(Schur_j(γ)) = ", norm_diff/opnorm(Schur_j))
        norm_diff = opnorm(Schur_m-Schur_jj)
        println("Difference Schur_m(γ) - Schur_jj(γ) = ", norm_diff)
        println("Relative difference opnorm(Schur_m(γ) - Schur_jj(γ))/opnorm(Schur_jj(γ)) = ", norm_diff/opnorm(Schur_jj))
    end
    println("\n--- End Schur-complement of native-WEP against MATLAB ---\n")
end

###########################################################
# Test the FFT-based Sylvester solver against naive solver
function fft_debug_mateq(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging FFT-Sylvester ---\n")
    γ = -rand(ComplexF64)
    gamma = γ
    C = rand(ComplexF64, nz, nx);
    waveguide = "JARLEBRING"


        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end
        println("  -- Matlab printouts start --")
        WEP_path = string(@__DIR__, "/WEP-ref")
         mat"""
            addpath($WEP_path)
            nxx = double($nx);
            nzz = double($nz);
            options = [];
            options.delta = double($delta);
            options.wg = $waveguide_str;
            nep = nep_wg_generator(nxx, nzz, options);

            CC = double($C);
            sigma = double($gamma);

            kk = mean(nep.K(:));
            K = nep.K - kk;

            $X_m = fft_wg( CC, sigma, kk, nep.hx, nep.hz );
         """
        println("  -- Matlab printouts end --")


    K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
    Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)

    k_bar = mean(K)

    A = Matrix(Dzz + 2*γ*Dz + (γ^2+k_bar)*I)
    B = complex(Matrix(Dxx))

    println("\nBuilt-in Sylvester solver (X_jj)")
    X_jj = @time sylvester(A,B,-C)
    println("Relative residual norm = ", opnorm(A*X_jj+X_jj*B-C)/opnorm(C))

    println("FFT-based Sylvester solver for WG (X_j)")
    X_j::Matrix{ComplexF64} = copy(C)
    @time solve_wg_sylvester_fft!( X_j, γ, k_bar, hx, hz )
    println("Relative residual norm = ", opnorm(A*X_j+X_j*B-C)/opnorm(C))

    println("MATLAB implemented FFT-based Sylvester solver for WG")
    println("Relative residual norm = ", opnorm(A*X_m+X_m*B-C)/opnorm(C))

    println("\nRelative difference opnorm(X_m - X_j)/opnorm(X_j) = ", opnorm(X_m - X_j)/opnorm(X_j))
    println("Relative difference opnorm(X_m - X_jj)/opnorm(X_jj) = ", opnorm(X_m - X_jj)/opnorm(X_jj))
    println("Relative difference opnorm(X_j - X_jj)/opnorm(X_j) = ", opnorm(X_j - X_jj)/opnorm(X_j))

    println("\n--- End FFT-Sylvester ---\n")
end


###########################################################
# Test the Sylvester SMW
function debug_Sylvester_SMW_WEP(nx::Integer, nz::Integer, delta::Number, N::Integer)
    println("\n\n--- Debugging Sylvester SMW for WEP ---\n")

    γ = -rand() - 1im*rand()
    gamma = γ
    σ = γ

    C = rand(Float64, nz, nx)

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing Sylvester SMW for waveguide: ", waveguide)

        nep = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "wEp", delta = delta)

        println("  Generate SMW-matrix")
        M_j = @time generate_smw_matrix(nep, N, σ)
        M_jj = Matrix(M_j)

        println("  Solve SMW-problem")
        X_j = @time solve_smw(nep, M_j, C, σ)

        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = string(@__DIR__, "/WEP-ref")
        mat"""
            addpath($WEP_path)
            nxx = double($nx);
            nzz = double($nz);
            nn = nzz;
            NN = double($N);
            CC = double($C);
            option = [];
            options.delta = double($delta);
            options.wg = $waveguide_str;
            nep = nep_wg_generator(nxx, nzz, options);

            sigma = double($gamma);

            %OBS: The below (MATLAB) code is copied from 'main.m' in the public code
            % CONSTRUCT THE PRECONDITIONER
            kk = mean(nep.K(:));
            K = nep.K - kk;


            Linv=@(X) fft_wg( X, sigma, kk, nep.hx, nep.hz );
            Pm_inv=@(x) -nep.Pm_inv(sigma, x);              %OBS! The minus sign!
            Pp_inv=@(x) -nep.Pp_inv(sigma, x);              %OBS! The minus sign!

            dd1 = nep.d1/nep.hx^2;
            dd2 = nep.d2/nep.hx^2;
            $M_m = generate_smw_matrix( nn, NN, Linv, dd1, dd2, Pm_inv, Pp_inv, K, false );

            $X_m = solve_smw( $M_m, CC, Linv, dd1, dd2, Pm_inv, Pp_inv, K );
        """
        println("  -- Matlab printouts end --")

        println("    Difference SMW_matrix_m - SMW_matrix_j = ", opnorm(M_m - M_jj))
        println("    Relative difference opnorm(SMW_matrix_m - SMW_matrix_j)/opnorm(SMW_matrix_j) = ", opnorm(M_m - M_jj)/opnorm(M_jj))
        println("\n    X is the solution to a Sylvester-SMW system")
        println("    Difference X_m - X_j = ", opnorm(X_m - X_j))
        println("    Relative difference opnorm(X_m - X_j)/opnorm(X_j) = ", opnorm(X_m - X_j)/opnorm(X_j))




    end
    println("\n--- End Debugging Sylvester SMW for WEP ---\n")
end


###########################################################
# Test convergence for the preconditioner
function debug_WEP_FD_preconditioner(delta::Number)
    println("\n\n--- Debugging WEP_FD preconditioner ---\n")

    nz = 3*3*5
    nx = nz + 4

    γ = -rand(ComplexF64)
    gamma = γ


    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Preconditioner in GMRES for waveguide: ", waveguide)

        nep = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "WeP", delta = delta)


        bb = rand(ComplexF64, nx*nz)
        precond = wep_generate_preconditioner(nep, 45, γ)
        Schur_fun = SchurMatVec(nep, γ)
        bbb = ldiv!(precond, copy(Schur_fun*bb))
        println("    Preconditioner * (Schur complement * b): Relative residual norm = ", norm(bbb - bb)/norm(bb))
        bb = rand(ComplexF64, nx*nz)
        bbb = ldiv!(precond, copy(bb))
        bbb = Schur_fun * bbb
        println("    Schur complement * (Preconditioner * b): Relative residual norm = ", norm(bbb - bb)/norm(bb))


        b = rand(ComplexF64, nx*nz+2*nz)
        for N = [1, 3, 9, 15, 45]
            println("    Testing for n = ", nz, " and N = ", N)
            precond = @time wep_generate_preconditioner(nep, N, γ)

            gmres_kwargs = ((:maxiter,100), (:restart,100), (:log,true), (:Pl,precond), (:reltol, 1e-13), (:verbose,true))
            linsolvercreator=GalleryWaveguide.WEPLinSolverCreator(solver_type=:gmres,kwargs=gmres_kwargs)
            wep_solver = create_linsolver(linsolvercreator, nep, γ)
            x = lin_solve(wep_solver, b)

            println("    Relative residual norm = ", norm(compute_Mlincomb(nep, gamma, x) - b)/norm(b))

        end


    end
    println("\n--- End WEP_FD preconditioner  ---\n")
end


###########################################################
# Test to compute eigenvalues and eigenvectors. Test preconditioned RESINV as in Ringh for native-WEP and MATLAB code, and pre-factorized-matrix RESINV for SPMF.
function matlab_debug_eigval_comp_WEP_FD_and_SPMF(nz::Integer, N::Integer, delta::Number)
    println("\n\n--- Debugging eigenvalue computations for WEP_FD and SPMF against MATLAB ---\n")

    nx = nz + 4;

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing eigenvalue computations for waveguide: ", waveguide)

        nep_j_WEPFD =    nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "WeP", delta = delta)
        nep_j_SPMF  =    nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "SpmF", delta = delta)
        nep_j_SPMF_pre = nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = waveguide, neptype = "SpmF_prE", delta = delta)

        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
            γ = -3-3.5im
        else
            waveguide_str = waveguide
            γ = -0.015-5.1im
        end
        gamma = γ


        println("    Generating preconditioner")
        precond = @time wep_generate_preconditioner(nep_j_WEPFD, N, γ)
        gmres_kwargs = ((:maxiter,100), (:restart,100), (:log,false), (:Pl,precond), (:tol, 1e-13))

        println("    Compute for WEP_FD")
        eigval_j_WEPFD = NaN
        eigvec_j_WEPFD = NaN
        try
            eigval_j_WEPFD, eigvec_j_WEPFD = @time resinv(nep_j_WEPFD, logger=1, λ=γ, maxit = 30,
                                                          tol = 1e-10, v=ones(ComplexF64,nx*nz+2*nz), c=zeros(ComplexF64,nx*nz+2*nz),
                                                          linsolvercreator=GalleryWaveguide.WEPLinSolverCreator(solver_type=:gmres,kwargs=gmres_kwargs)
                                                          );
        catch err
            # Only catch NoConvergence
            isa(err, NoConvergenceException) || rethrow(err)
            println("No convergence because: ", err.msg)
            # still access the approximations
            eigval_j_WEPFD = err.λ
            eigvec_j_WEPFD = err.v
        end


        println("    Compute for SPMF")
        eigval_j_SPMF = NaN
        eigvec_j_SPMF = NaN
        try
            eigval_j_SPMF, eigvec_j_SPMF = @time resinv(nep_j_SPMF, logger=1, λ=γ, maxit = 30, tol = 1e-10, v=ones(ComplexF64,nx*nz+2*nz), c=zeros(ComplexF64,nx*nz+2*nz))
        catch err
            # Only catch NoConvergence
            isa(err, NoConvergenceException) || rethrow(err)
            println("No convergence because: ", err.msg)
            # still access the approximations
            eigval_j_SPMF = err.λ
            eigvec_j_SPMF = err.v
        end

        println("    Compute for SPMF (pre-Schur-fact)")
        eigval_j_SPMF_pre = NaN
        eigvec_j_SPMF_pre = NaN
        try
            eigval_j_SPMF_pre, eigvec_j_SPMF_pre = @time resinv(nep_j_SPMF_pre, logger=1, λ=γ, maxit = 30, tol = 1e-10, v=ones(ComplexF64,nx*nz+2*nz), c=zeros(ComplexF64,nx*nz+2*nz))
        catch err
            # Only catch NoConvergence
            isa(err, NoConvergenceException) || rethrow(err)
            println("No convergence because: ", err.msg)
            # still access the approximations
            eigval_j_SPMF_pre = err.λ
            eigvec_j_SPMF_pre = err.v
        end


        println("  -- Matlab printouts start --")
        WEP_path = string(@__DIR__, "/WEP-ref")
        mat"""
            addpath($WEP_path)

            nzz = double($nz);
            lambda = double($gamma);
            NN = double($N);

            [$eigval_m, $eigvec_m] = main_func_WEP(nzz, NN, lambda, $delta, $waveguide_str)

        """
        println("  -- Matlab printouts end --")

        println("    Difference between WEP_FD and SPMF computed eigenvalue = ", abs(eigval_j_WEPFD - eigval_j_SPMF))
        println("    Relative difference between WEP_FD and SPMF computed eigenvalue = ", abs(eigval_j_WEPFD - eigval_j_SPMF)/abs(eigval_j_WEPFD))
        println("    Difference between WEP_FD and SPMF computed eigenvectors = ", norm(eigvec_j_WEPFD/eigvec_j_WEPFD[1] - eigvec_j_SPMF/eigvec_j_SPMF[1]))
        println("    Relative difference between WEP_FD and SPMF computed eigenvalue = ", norm(eigvec_j_WEPFD/eigvec_j_WEPFD[1] - eigvec_j_SPMF/eigvec_j_SPMF[1])/norm(eigvec_j_WEPFD/eigvec_j_WEPFD[1]))
        println("")
        println("    Difference between WEP_FD and MATLAB computed eigenvalue = ", abs(eigval_j_WEPFD - eigval_m))
        println("    Relative difference between WEP_FD and MATLAB computed eigenvalue = ", abs(eigval_j_WEPFD - eigval_m)/abs(eigval_m))
        println("    Difference between WEP_FD and MATLAB computed eigenvectors = ", norm(eigvec_j_WEPFD/eigvec_j_WEPFD[1] - eigvec_m/eigvec_m[1]))
        println("    Relative difference between WEP_FD and MATLAB computed eigenvalue = ", norm(eigvec_j_WEPFD/eigvec_j_WEPFD[1] - eigvec_m/eigvec_m[1])/norm(eigvec_m/eigvec_m[1]))
        println("")
        println("    Difference between MATLAB and SPMF computed eigenvalue = ", abs(eigval_m - eigval_j_SPMF))
        println("    Relative difference between MATLAB and SPMF computed eigenvalue = ", abs(eigval_m - eigval_j_SPMF)/abs(eigval_m))
        println("    Difference between MATLAB and SPMF computed eigenvectors = ", norm(eigvec_m/eigvec_m[1] - eigvec_j_SPMF/eigvec_j_SPMF[1]))
        println("    Relative difference between MATLAB and SPMF computed eigenvalue = ", norm(eigvec_m/eigvec_m[1] - eigvec_j_SPMF/eigvec_j_SPMF[1])/norm(eigvec_m/eigvec_m[1]))
        println("")
        println("    Difference between SPMF and SPMF-pre computed eigenvalue = ", abs(eigval_j_SPMF_pre - eigval_j_SPMF))
        println("    Relative difference between SPMF and SPMF-pre computed eigenvalue = ", abs(eigval_j_SPMF_pre - eigval_j_SPMF)/abs(eigval_j_SPMF))
        println("    Difference between SPMF and SPMF-pre computed eigenvectors = ", norm(eigvec_j_SPMF_pre/eigvec_j_SPMF_pre[1] - eigvec_j_SPMF/eigvec_j_SPMF[1]))
        println("    Relative difference between SPMF and SPMF-pre computed eigenvalue = ", norm(eigvec_j_SPMF_pre/eigvec_j_SPMF_pre[1] - eigvec_j_SPMF/eigvec_j_SPMF[1])/norm(eigvec_j_SPMF/eigvec_j_SPMF[1]))
        println("")
        println("    WEP_FD converged to the eigenvalue = ", eigval_j_WEPFD)

    end
    println("\n--- End eigenvalue computations of WEP_FD and SPMF against MATLAB ---\n")
end



end
