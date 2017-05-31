
using MATLAB
using Gallery

#OBS: Need 'using Waveguide' in order to run the debug tests
# Gallery is not exporting this module

export matlab_debug_WEP_FD #ONLY FOR DEBUGGING
export matlab_debug_full_matrix_WEP_FD_SPMF #ONLY FOR DEBUGGING
export debug_sqrtm_schur #ONLY FOR DEBUGGING
export fft_debug_mateq #ONLY FOR DEBUGGING
export debug_sqrt_derivative #ONLY FOR DEBUGGING
export debug_Mlincomb_FD_WEP #ONLY FOR DEBUGGING
export debug_Sylvester_SMW_WEP #ONLY FOR DEBUGGING



########### SOME REFERENCE IMPLEMENTATIONS ################
###########################################################

# Generate P-matrix
# Is the lower right part of the system matrix, from the DtN maps Jarlebring-(1.5)(1.6) and Ringh-(2.4)(2.8)
function generate_P_matrix(nz::Integer, hx, Km, Kp)

    R, Rinv = generate_R_matrix(nz::Integer)
    const p = (nz-1)/2;

    # Constants from the problem
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


######################## DEBUG ############################
###########################################################

# Test the basic generated FD matrices against MATLAB code
function matlab_debug_WEP_FD(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging Matrices FD against MATLAB ---\n")
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

        K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        P, p_P = generate_P_matrix(nz, hx, Km, Kp)
        P_j = zeros(Complex128, 2*nz,2*nz)
        Iz = eye(2*nz, 2*nz)
        for i = 1:2*nz
            P_j[:,i] = P(γ, Iz[:,i])
        end

        R, Rinv = generate_R_matrix(nz)
        S = generate_S_function(nz, hx, Km, Kp)
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
        WEP_path = pwd() * "/../matlab/WEP"
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
    println("\n--- End Matrices FD against MATLAB ---\n")
end


###########################################################
# Test the full generated system-matrix against against MATLAB code
function matlab_debug_full_matrix_WEP_FD_SPMF(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging Full Matrix FD SPMF against MATLAB ---\n")
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
        WEP_path = pwd() * "/../matlab/WEP"
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
    println("\n--- End Full Matrix FD SPMF against MATLAB ---\n")
end


###########################################################
# Test the square-root implementations. Both reference and negative-imaginary-part.
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
# Test computation of derivative <d> of sqrt(ax^2 + bx + c) in <x>
function debug_sqrt_derivative()
    println("\n\n--- Debugging derivatives of square root of polynomials against MATLAB symbolic ---\n")
    a = 2*rand()
    b = 2*pi*rand(Complex128)
    c = 1.67*rand()
    d_vec = [0 1 2 3 4 11 19 20 21 22 30 35 45 60] #Factorial for Int64 overflows at 21!
    x = 25*rand()

        WEP_path = pwd() * "/../matlab/WEP"
        println("  -- Matlab printouts start --")
        @mput a b c d_vec x WEP_path
        @matlab begin
            addpath(WEP_path)
            der_val = sqrt_derivative_test(a,b,c, d_vec, x);
        @matlab end
        @mget der_val
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
        warn("This debug is 'naive' and might be slow for the discretization used.")
    end


    γ = -rand() - 1im*rand()
    gamma = γ
    n = nx*nz+2*nz
    V = rand(Complex128, n, 4)

    for waveguide = ["TAUSCH", "JARLEBRING"]
        println("\n")
        println("Testing full matrix M for waveguide: ", waveguide)

        nep_j = nep_gallery("waveguide", nx, nz, waveguide, "fD", "weP", delta)
        nep_SPMF = nep_gallery("waveguide", nx, nz, waveguide, "fD", "SpmF", delta)

        for d = [0 1 2 5]
            M_j = zeros(Complex128, n, n)
            for i = 1:n
                V = zeros(Complex128, n)
                V[i] = 1
                M_j[:,i] += compute_Mlincomb(nep_j, γ, V, [1], d)
            end

            M_SPMF = zeros(Complex128, n, n)
            for i = 1:n
                V = zeros(Complex128, n)
                V[i] = 1
                M_SPMF[:,i] += compute_Mlincomb(nep_SPMF, γ, V, [1], d)
            end
            println("Derivative d = ", d)
            println("    Difference M_SPMF(γ) - M(γ) = ", norm(full(M_SPMF-M_j)))
            println("    Relative difference norm(M_m(γ) - M(γ))/norm(M(γ)) = ", norm(full(M_SPMF-M_j))/norm(full(M_j)))

            v_j = compute_Mlincomb(nep_j, γ, V)
            v_SPMF = compute_Mlincomb(nep_SPMF, γ, V)
            println("    Relative difference on random set of 4 vectors = ", norm(v_j - v_SPMF)/norm(v_j))
        end


    end
    println("\n--- End Mlincomb (Full Matrix FD WEP-native-format against SPMF) ---\n")
end


###########################################################
# Test the FFT-based Sylvester solver against naive solver
function fft_debug_mateq(nx::Integer, nz::Integer, delta::Number)
    println("\n\n--- Debugging FFT-Sylvester ---\n")
    γ = -rand(Complex128)
    gamma = γ
    C = rand(Complex128, nz, nx);
    waveguide = "JARLEBRING"


        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end
        println("  -- Matlab printouts start --")
        WEP_path = pwd() * "/../matlab/WEP"
        @mput nx nz delta WEP_path waveguide_str gamma C
        @matlab begin
            addpath(WEP_path)
            nxx = double(nx)
            nzz = double(nz)
            options = struct
            options.delta = double(delta)
            options.wg = waveguide_str
            nep = nep_wg_generator(nxx, nzz, options)

            CC = double(C)
            sigma = double(gamma)

            kk = mean(nep.K(:));
            K = nep.K - kk;

            X_m = fft_wg( CC, sigma, kk, nep.hx, nep.hz );

        @matlab end
        @mget X_m
        println("  -- Matlab printouts end --")


    K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
    Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)

    k_bar = mean(K)

    A = full(Dzz + 2*γ*Dz + (γ^2+k_bar)*speye(Complex128, nz,nz))
    B = complex(full(Dxx))

    println("\nBuilt-in Sylvester solver (X_jj)")
    X_jj = @time sylvester(A,B,-C)
    println("Relative residual norm = ", norm(A*X_jj+X_jj*B-C)/norm(C))

    println("FFT-based Sylvester solver for WG (X_j)")
    X_j = @time solve_wg_sylvester_fft( C, γ, k_bar, hx, hz )
    println("Relative residual norm = ", norm(A*X_j+X_j*B-C)/norm(C))

    println("MATLAB implemented FFT-based Sylvester solver for WG")
    println("Relative residual norm = ", norm(A*X_m+X_m*B-C)/norm(C))

    println("\nRelative difference norm(X_m - X_j)/norm(X_j) = ", norm(X_m - X_j)/norm(X_j))
    println("Relative difference norm(X_m - X_jj)/norm(X_jj) = ", norm(X_m - X_jj)/norm(X_jj))
    println("Relative difference norm(X_j - X_jj)/norm(X_j) = ", norm(X_j - X_jj)/norm(X_j))

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

        nep = nep_gallery("waveguide", nx, nz, waveguide, "fD", "weP", delta)

        M_j = @time generate_smw_matrix(nep, N, σ)
        M_jj = full(M_j)

        X_j = @time solve_smw(nep, M_j, C, σ)

        if waveguide == "JARLEBRING"
            waveguide_str = "CHALLENGE"
        else
            waveguide_str = waveguide
        end

        println("  -- Matlab printouts start --")
        WEP_path = pwd() * "/../matlab/WEP"
        @mput nx nz N delta WEP_path waveguide_str gamma C
        @matlab begin
            addpath(WEP_path)
            nxx = double(nx)
            nzz = double(nz)
            nn = nz
            NN = double(N)
            CC = double(C)
            options = struct
            options.delta = double(delta)
            options.wg = waveguide_str
            nep = nep_wg_generator(nxx, nzz, options)

            sigma = double(gamma)

            #OBS: The below (MATLAB) code is copied from 'main.m' in the public code
            # CONSTRUCT THE PRECONDITIONER
            kk = mean(nep.K(:));
            K = nep.K - kk;


            eval("Linv=@(X) fft_wg( X, sigma, kk, nep.hx, nep.hz );")
            eval("Pm_inv=@(x) -nep.Pm_inv(sigma, x);              %OBS! The minus sign!")
            eval("Pp_inv=@(x) -nep.Pp_inv(sigma, x);              %OBS! The minus sign!")

            dd1 = nep.d1/nep.hx^2;
            dd2 = nep.d2/nep.hx^2;
            M_m = generate_smw_matrix( nn, NN, Linv, dd1, dd2, Pm_inv, Pp_inv, K, false );

            X_m = solve_smw( M_m, CC, Linv, dd1, dd2, Pm_inv, Pp_inv, K );

        @matlab end
        @mget M_m X_m
        println("  -- Matlab printouts end --")

        println("    Difference SMW_matrix_m - SMW_matrix_j = ", norm(M_m - M_jj))
        println("    Relative difference norm(SMW_matrix_m - SMW_matrix_j)/norm(SMW_matrix_j) = ", norm(M_m - M_jj)/norm(M_jj))
        println("\n    X is the solution to a Sylvester-SMW system")
        println("    Difference X_m - X_j = ", norm(X_m - X_j))
        println("    Relative difference norm(X_m - X_j)/norm(X_j) = ", norm(X_m - X_j)/norm(X_j))




    end
    println("\n--- End Debugging Sylvester SMW for WEP ---\n")
end
