
###########################################################
# Sylvester based preconditioner for the Waveguide eigenvalue problem
# Implementing the Preconditioner in the optimized way described in Ringh et al.

"""
    WEP_preconditioner
 A wrapper struct to be able to use the preconditioner. Update in GMRES for Julia v0.6 remove possibility to use a function directly.
"""
    struct WEP_preconditioner
        scratch_pad_for_FFT::Matrix{ComplexF64}
        scratch_pad_for_transpose::Matrix{ComplexF64}
        scratch_pad_for_Z::Matrix{ComplexF64}
        nep::WEP_FD
        M
        σ
        function WEP_preconditioner(nep::NEP, M, σ)
            nx = nep.nx
            nz = nep.nz
            A = zeros(ComplexF64, 2*(nx+1), nz)
            B = zeros(ComplexF64, nx, nz)
            C = zeros(ComplexF64, nz, nx)
            return new(A, B, C, nep, M, σ)
        end
    end


    function ldiv!(A::WEP_preconditioner, B)
        C = reshape(B, A.nep.nz, A.nep.nx)
        B[:] = vec( solve_smw( A.nep, A.M, C, A.σ, A.scratch_pad_for_FFT, A.scratch_pad_for_transpose, A.scratch_pad_for_Z) )
    end



# Ringh - Algorithm 2, step 10
"""
    wep_generate_preconditioner(nep::WEP_FD, N::Integer, σ)
 Given a nep of type WEP_FD, the number of domains in z-direction N, and a fixed shift σ,\\
 this computes a functor (struct acting like a function) that acts as a Preconditioner for the WEP.
"""
    function wep_generate_preconditioner(nep::WEP_FD, N::Integer, σ)
        M = generate_smw_matrix(nep, N, σ)
        return WEP_preconditioner(nep, M, σ)
    end


"""
    generate_smw_matrix( nep::WEP_FD, N::Integer, σ)
 Given a nep of type WEP_FD, this computes the Sylvester-SMW matrix with N domains in z-direction\\
 and for fixed shift σ.
"""
    function generate_smw_matrix(nep::WEP_FD, N::Integer, σ)
        # OBS: n = nz, and nz = nx + 4;
        if( (nep.nz+4) != nep.nx)
            error("This implementation requires nx = nz + 4. Provided NEP has nz = ", nep.nz, " and nx = ", nep.nx)
        end
        if( !isinteger(nep.nz/N) )
            error("This implementation is uniform in the blocking and therefore requires nz/N tobe an integer. Provided data is nz = ", nep.nz, " with N = ", N, " and hence nz/N = ", nep.nz/N, " which is not deemed to be numerically equal to an integer." )
        end

        nz::Integer = nep.nz
        nx::Integer = nep.nx

        dd1 = nep.d1/nep.hx^2;
        dd2 = nep.d2/nep.hx^2;

        # Sylvester solver
        scratch_pad_for_FFT::Matrix{ComplexF64} = zeros(ComplexF64, 2*(nx+1), nz)
        scratch_pad_for_transpose::Matrix{ComplexF64} = zeros(ComplexF64, nx, nz)
        scratch_pad_for_Z::Matrix{ComplexF64} = zeros(ComplexF64, nz, nx)
        Linv! = function(rhs)
            return solve_wg_sylvester_fft!(rhs, σ, nep.k_bar, nep.hx, nep.hz, scratch_pad_for_FFT, scratch_pad_for_transpose, scratch_pad_for_Z)
        end

        # OBS: MINUS sign as in Ringh - (4.10)
        Pm(v) = -P_inv_m(nep, σ, v)
        Pp(v) = -P_inv_p(nep, σ, v)

        return generate_smw_matrix(nz, N, Linv!, dd1, dd2, Pm, Pp, nep.K )
    end


"""
    solve_smw(nep::WEP_FD, M, C, σ, scratch_pad_for_FFT, scratch_pad_for_transpose, scratch_pad_for_Z)
 Given a nep of type WEP_FD, an SMW-system matrix M computed with shift σ, and a right hand side C,\\
 This computes the solution to the SMW-matrix equation.
"""
    solve_smw(nep::WEP_FD, M, C, σ) = solve_smw(nep::WEP_FD, M, C, σ, zeros(ComplexF64, 2*(size(C,2)+1), size(C,1)), zeros(ComplexF64, size(C,2), size(C,1)), zeros(ComplexF64, size(C,1), size(C,2)))

    function solve_smw(nep::WEP_FD, M, C, σ, scratch_pad_for_FFT, scratch_pad_for_transpose, scratch_pad_for_Z)

        C_copy::Matrix{ComplexF64} = copy(C) #Make sure it is complex since that is how FFT works. Also take copy since WG_FFT solver works in place.

        nz::Integer = nep.nz
        nx::Integer = nep.nx

        dd1 = nep.d1/nep.hx^2;
        dd2 = nep.d2/nep.hx^2;

        # Sylvester solver
        Linv! = function(rhs)
            return solve_wg_sylvester_fft!(rhs, σ, nep.k_bar, nep.hx, nep.hz, scratch_pad_for_FFT, scratch_pad_for_transpose, scratch_pad_for_Z)
        end

        # OBS: MINUS sign as in Ringh - (4.10)
        Pm(v) = -P_inv_m(nep, σ, v)
        Pp(v) = -P_inv_p(nep, σ, v)

        return solve_smw(M, C_copy, Linv!, dd1, dd2, Pm, Pp, nep.K)
    end


#FFT diagonalization and solution to WEP matrix equation.
#Ringh - Section 5.3
"""
    solve_wg_sylvester_fft( C, λ, k_bar, hx, hz, scratch_pad_for_FFT, scratch_pad_for_transpose )
 Solves the Sylvester equation for the WEP, with C as right hand side.
 OBS: Writes the solution to the varaible C.
 Last three arguments (scratch pad:s) are optional and there to allow reuse of memory allocation. They will be overwritten in the process!
"""
solve_wg_sylvester_fft!( C, λ, k_bar, hx, hz) = solve_wg_sylvester_fft!( C, λ, k_bar, hx, hz, zeros(ComplexF64, 2*(size(C,2)+1), size(C,1)), zeros(ComplexF64, size(C,2), size(C,1)), zeros(ComplexF64, size(C,1), size(C,2)))

function solve_wg_sylvester_fft!( C, λ, k_bar, hx, hz, scratch_pad_for_FFT, scratch_pad_for_transpose, scratch_pad_for_Z)

    nz = size(C,1)
    nx = size(C,2)

    alpha = λ^2+k_bar;

    #eigenvalues of A
    v=zeros(ComplexF64, nz);   v[1]=-2;    v[2]=1;     v[nz]=1;   v=v/(hz^2);
    w=zeros(ComplexF64, nz);   w[2]=1;     w[nz]=-1;   w=w*(λ/hz);
    D = fft(v+w) .+ alpha

    # eigenvalues of B = Dxx
    S = -((4+0.0im)/hx^2) * sin.(pi*(1:nx)/(2*(nx+1))).^2
#    S=S.'


    # change variables
    adjoint!(scratch_pad_for_transpose, C)  # scratch_pad_for_transpose = C'
    adjoint!(C, Wh(scratch_pad_for_transpose, scratch_pad_for_FFT )) # C = (Wh(C'))'
    Vh!( C )  #In effect: C = Vh( Wh(C')' )

    # solve the diagonal matrix equation
    scratch_pad_for_Z[:,:] .= 0.0im
    for k=1:nx
        scratch_pad_for_Z[:,k] += C[:,k] ./ (D .+ S[k])
    end


    # change variables
    adjoint!(scratch_pad_for_transpose, scratch_pad_for_Z)  # scratch_pad_for_transpose = Z'
    adjoint!(C, W(scratch_pad_for_transpose, scratch_pad_for_FFT ))  #C = (W(Z'))'
    V!( C ) #In effect: C = V( W(Z')' )

    return C

end


# Start: Auxiliary computations of eigenvector actions using FFT
#Connects to the solving of the Sylvester equation by FFT-diagonalization
#Ringh - Section 5.3
    function V!(X::Matrix{ComplexF64})
    # Compute the action of the eigenvectors of A = Dzz + Dz + c*I
        nx::ComplexF64 = size(X,2)
        fft!(X,1)#/sqrt(nx)
        rmul!(X, (1+0.0im)/sqrt(nx))
    end

    function Vh!(X::Matrix{ComplexF64})
    # Compute the action of the transpose of the eigenvectors of A = Dzz + Dz + c*I
        nx::ComplexF64 = size(X,2)
        ifft!(X,1)#* sqrt(nx)
        rmul!(X, (1+0.0im)*sqrt(nx))
    end

    function W(X::Matrix{ComplexF64}, scratch_pad::Matrix{ComplexF64})
    #W Compute the action of the matrix W
    #   W is the matrix of the eigenvectors of the second derivative Dxx
    #   W*X can be computed with FFTs

        nz::ComplexF64 = size(X,1)
        WX::Matrix{ComplexF64} = F(X,scratch_pad)-Fh(X,scratch_pad)

        rmul!(WX, (1.0im/2.0) * 1/sqrt((nz+1)/2.0) )
        return WX

    end

    function Wh(X::Matrix{ComplexF64}, scratch_pad::Matrix{ComplexF64})
    #Wh Compute the action of the matrix Wh
    #   Wh is the transpose of the matrix of the eigenvectors of the second derivative Dxx
    #   Wh*X can be computed with FFTs

        nz::ComplexF64 = size(X,1)
        WX::Matrix{ComplexF64} = F(X,scratch_pad)-Fh(X,scratch_pad)

        rmul!(WX, (1.0im/2.0) * 1/sqrt((nz+1)/2.0) )
        return WX

    end

    function F(v::Matrix{ComplexF64}, scratch_pad::Matrix{ComplexF64})
    #F is an auxiliary function for W and Wh

        m=size(v,2)
        n=size(v,1) + 1

        scratch_pad[:,:] .= 0.0im
        scratch_pad[2:n,:] = v
        fft!(scratch_pad,1)
        return scratch_pad[2:n,:]
    end

    function Fh( v::Matrix{ComplexF64}, scratch_pad::Matrix{ComplexF64})
    #Fh is an auxiliary function for W and Wh

        m=size(v,2)
        n=size(v,1) + 1

        scratch_pad[:,:] .= 0.0im
        scratch_pad[2:n,:] = v
        ifft!(scratch_pad,1)
        return scratch_pad[2:n,:]*2*n

    end
# End: Auxiliary computations of eigenvector actions using FFT


###########################################################
# Sylvester SMW
# Ringh - Section 4
"""
    generate_smw_matrix(n::Integer, N::Integer, Linv!::Function, dd1, dd2, Pm, Pp, K)
 Computes the SMW matrix for the Sylvester SMW on rectangular domains, with n points in z-direction, n+4 points in x-direction,\\
 N domains in z-direction, N+4 domains in x-direction, and fixed shift.\\
 Obs: Linv! works in place on the matrix rhs
"""
function generate_smw_matrix(n::Integer, N::Integer, Linv!::Function, dd1, dd2, Pm::Function, Pp::Function, K::Union{Matrix{ComplexF64}, Matrix{Float64}})

    # OBS: n = nz, and nz = nx + 4
    nz::Integer = n
    nx::Integer = n + 4


    L::Integer = n/N             # Number of points in one dimanesion of the regions
    LL::Integer = L*L            # Number of points in the "square interior regions"
    mm::Integer = (N^2 + 4*N)    # Number of elements in SMW-matrix

    # block extract index
    II = function (i::Integer) # z-direction, always equal
        return (i-1)*L+1:i*L
    end
    JJ = function(j::Integer) #x-direction, different if boundary or interior
        return ((j-3)*L+1:(j-2)*L) .+ 2
    end
    JJ_2 = function(j::Integer)
        return (j==1)*1 + (j==2)*2 + (j==N+3)*(n+3) + (j==N+4)*(n+4)
    end

    # convert a single index k to two indeces
    # reading the matrix X by rows; left to right and top to down
    k2ij = function(k::Integer)
        j::Integer = rem(k,N+4)+(rem(k,N+4)==0)*(N+4)
        i::Integer = (k-j)/(N+4)+1
        return(i,j)
    end

    # compute matrix M
    M::Matrix{ComplexF64} = zeros(ComplexF64, mm, mm)

    EEk::Matrix{ComplexF64} = zeros(ComplexF64, nz, nx)
    ek::Vector{ComplexF64} = zeros(ComplexF64, nz)

    for k=1:mm

        i,j = k2ij(k)

        # Ek tilde
        EEk[:,:] .= 0.0im
        if j==1
            EEk[II(i), JJ_2(j)] = K[II(i), JJ_2(j)]
            ek[:] .= 0.0im
            ek[II(i)] .= dd1
            EEk[:, 1] += Pm(ek)
        elseif j==2
            EEk[II(i), JJ_2(j)] = K[II(i), JJ_2(j)]
            ek[:] .= 0.0im
            ek[II(i)] .= dd2
            EEk[:, 1] += Pm(ek)
        elseif j==N+4
            EEk[II(i),JJ_2(j)] = K[II(i),JJ_2(j)]
            ek[:] .= 0.0im
            ek[II(i)] .= dd1
            EEk[:, nx] += Pp(ek)
        elseif j==N+3
            EEk[II(i), JJ_2(j)] = K[II(i), JJ_2(j)]
            ek[:] .= 0.0im
            ek[II(i)] .= dd2
            EEk[:, nx] += Pp(ek)
        else
            EEk[II(i), JJ(j)] = K[II(i), JJ(j)]
        end

        # Sylvester solve of E tilde
        Linv!(EEk) # In effect EEk = Linv!(EEk) = Fk

        # Build this matrix element

        for kk=1:mm
            i,j = k2ij(kk)
            # evaluate the linear functional
            if((j==1)||(j==2)||(j==N+3)||(j==N+4))
                M[kk, k] = sum(sum(EEk[II(i), JJ_2(j)] ))/L
            else
                M[kk, k] = sum(sum(EEk[II(i), JJ(j)] ))/LL
            end
        end
    end

    M += Matrix{ComplexF64}(I, mm, mm)

    return factorize(M)

end



"""
    solve_smw( M, C, Linv!::Function, dd1, dd2, Pm, Pp, K)
 Solves the matrix equation SMW system and computes the solution, on rectangular domains.\n
 With SMW-system matrix M, right hand side C, and matrix equation solver Linv! which was used to compute M.
 Obs: Linv! works in place on the matrix rhs
"""
function solve_smw( M, C::Matrix{ComplexF64}, Linv!::Function, dd1, dd2, Pm::Function, Pp::Function, K::Union{Matrix{ComplexF64}, Matrix{Float64}})

    mm::Integer = size(M,1)
    N::Integer = sqrt(mm+4)-2      #OBS: N^2 + 4N = length(M)

    nz::Integer = size(C,1)
    nx::Integer = size(C,2)
    n::Integer = nz
    L::Integer = nz/N
    LL::Integer = L*L

    # block extract index
    II = function (i::Integer) # z-direction, always equal
        return (i-1)*L+1:i*L
    end
    JJ = function(j::Integer) #x-direction, different if boundary or interior
        return ((j-3)*L+1:(j-2)*L) .+ 2
    end
    JJ_2 = function(j::Integer)
        return (j==1)*1 + (j==2)*2 + (j==N+3)*(n+3) + (j==N+4)*(n+4)
    end

    # convert a single index k to two indeces
    # reading the matrix X by rows; left to right and top to down
    k2ij = function(k::Integer)
        j::Integer = rem(k,N+4)+(rem(k,N+4)==0)*(N+4)
        i::Integer = (k-j)/(N+4)+1
        return(i,j)
    end

    # compute the right hand side
    Linv!(C) # Now C = LinvC

    b=zeros(ComplexF64,mm)
    for k=1:mm
        i, j = k2ij(k)
        # evaluate the linear functional
        if((j==1)||(j==2)||(j==N+3)||(j==N+4))
            b[k] = sum(sum( C[II(i), JJ_2(j)] ))/L
        else
            b[k] = sum(sum( C[II(i), JJ(j)] ))/LL
        end
    end

    # solve for the coefficients
    alpha = M\b;

    # build the solution
    Y::Matrix{ComplexF64} = zeros(ComplexF64, nz, nx);
    ek::Vector{Float64} = zeros(ComplexF64, nz);
    for k=1:mm

        i, j = k2ij(k)

        # Ek tilde implicitly created
        if j==1
            Y[II(i), 1] += alpha[k]*K[II(i), 1]
            ek[:] .= 0.0im
            ek[II(i)] .= dd1
            Y[:, 1] += alpha[k]*Pm(ek)
        elseif j==2
            Y[II(i), 2] += Y[II(i),2] + alpha[k]*K[II(i), 2]
            ek[:] .= 0.0im
            ek[II(i)] .= dd2
            Y[:, 1] += alpha[k]*Pm(ek)
        elseif j==N+4
            Y[II(i), nx] += alpha[k]*K[II(i), nx]
            ek[:] .= 0.0im
            ek[II(i)] .= dd1
            Y[:, nx] += alpha[k]*Pp(ek)
        elseif j==N+3
            Y[II(i), nx-1] += alpha[k]*K[II(i), nx-1]
            ek[:] .= 0.0im
            ek[II(i)] .= dd2
            Y[:, nx] += alpha[k]*Pp(ek)
        else
            Y[II(i),JJ(j)] += alpha[k]*K[II(i), JJ(j)]
        end

    end
    Linv!(Y) # Now Y = LinvY

    return C-Y # Which is LinvC - LinvY

end
