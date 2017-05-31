
###########################################################
# Sylvester based preconditioner for the Waveguide eigenvalue problem
# Implementing the Preconditioner in the optimized way described in Ringh et al.

import IterativeSolvers.solve
export solve

#FFT diagonalization and solution to WEP matrix equation.
#Ringh - Section 5.3
    """
    solve_wg_sylvester_fft( C, λ, k_bar, hx, hz )
 Solves the Sylvester equation for the WEP, with C as right hand side.
"""
function solve_wg_sylvester_fft( C, λ, k_bar, hx, hz )

    nz = size(C,1)
    nx = size(C,2)

    alpha = λ^2+k_bar;

    #eigenvalues of A
    v=zeros(nz);   v[1]=-2;    v[2]=1;     v[nz]=1;   v=v/(hz^2);
    w=zeros(nz);   w[2]=1;     w[nz]=-1;   w=w*(λ/hz);
    D=fft(v+w)+alpha;

    # eigenvalues of B = Dxx
    S = -(4/hx^2) * sin(pi*(1:nx)/(2*(nx+1))).^2
#    S=S.'

    # solve the diagonal matrix equation
    Z=zeros(Complex128,nz,nx)

    CC = Vh!( Wh(C' )' )

    for k=1:nx
        Z[:,k] += CC[:,k]./(D+S[k])
    end


    # change variables
    return V!((W(Z'))')

end


# Start: Auxiliary computations of eigenvector actions using FFT
#Connects to the solving of the Sylvester equation by FFT-diagonalization
#Ringh - Section 5.3
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

        WX = (1.0im/2.0)*(F(X)-Fh(X))

        nz = size(X,1)
        return WX/sqrt((nz+1)/2.0)

    end

    function Wh( X )
    #Wh Compute the action of the matrix Wh
    #   Wh is the transpose of the matrix of the eigenvectors of the second derivative Dxx
    #   Wh*X can be computed with FFTs

        WX = (Fh(X)-F(X))/(2.0im)

        nz = size(X,1)
        return WX/sqrt((nz+1)/2.0)

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


###########################################################
# Sylvester SMW
# Ringh - Section 4
    """
    generate_smw_matrix(n::Integer, N::Integer, Linv::Function, dd1, dd2, Pm, Pp, K)
 Computes the SMW matrix for the Sylvester SMW on rectangular domains, with n points in z-direction, n+4 points in x-direction,\\
 N domains in z-direction, N+4 domains in x-direction.
 and fixed shift.
"""
function generate_smw_matrix(n::Integer, N::Integer, Linv::Function, dd1, dd2, Pm::Function, Pp::Function, K)

    # OBS: n = nz, and nz = nx + 4
    const nz::Integer = n
    const nx::Integer = n + 4

        
    const L::Integer = n/N             # Number of points in one dimanesion of the regions
    const LL::Integer = L*L            # Number of points in the "square interior regions"
    const mm::Integer = (N^2 + 4*N)    # Number of elements in SMW-matrix

    # block extract index
    II = function (i::Integer) # z-direction, always equal
        (i-1)*L+1:i*L
    end
    JJ = function(j::Integer) #x-direction, different if boundary or interior
        ((j-3)*L+1:(j-2)*L) + 2
    end
    JJ_2 = function(j::Integer)
        (j==1)*1 + (j==2)*2 + (j==N+3)*(n+3) + (j==N+4)*(n+4)
    end

    # convert a single index k to two indeces
    # reading the matrix X by rows; left to right and top to down
    k2ij = function(k::Integer)
        j::Integer = rem(k,N+4)+(rem(k,N+4)==0)*(N+4)
        i::Integer = (k-j)/(N+4)+1
        return(i,j)
    end

    # compute matrix M
    M::Array{Complex128,2} = zeros(Complex128, mm, mm)

    EEk::Array{Complex128,2} = zeros(Complex128, nz, nx)
    ek::Array{Complex128,1} = zeros(Complex128, nz)

    for k=1:mm

        i,j = k2ij(k);

        # Ek tilde
        EEk = 0*EEk;
        if (j==1)
            EEk[II(i), JJ_2(j)] = K[II(i), JJ_2(j)]
            ek = 0*ek
            ek[II(i)] = dd1
            EEk[:, 1] += Pm(ek)
        elseif (j==2)
            EEk[II(i), JJ_2(j)] = K[II(i), JJ_2(j)]
            ek = 0*ek
            ek[II(i)] = dd2
            EEk[:, 1] += Pm(ek)
        elseif (j==N+4)
            EEk[II(i),JJ_2(j)] = K[II(i),JJ_2(j)]
            ek = 0*ek
            ek[II(i)] = dd1
            EEk[:, nx] += Pp(ek)
        elseif (j==N+3)
            EEk[II(i), JJ_2(j)] = K[II(i), JJ_2(j)]
            ek = 0*ek
            ek[II(i)] = dd2
            EEk[:, nx] += Pp(ek)
        else
            EEk[II(i), JJ(j)] = K[II(i), JJ(j)]
        end

        # Sylvester solve of E tilde
        Fk = Linv(EEk);

        # Build this matrix element

        for kk=1:mm
            i,j = k2ij(kk);
            # evaluate the linear functional
            if((j==1)||(j==2)||(j==N+3)||(j==N+4))
                M[kk, k] = sum(sum(Fk[II(i), JJ_2(j)] ))/L
            else
                M[kk, k] = sum(sum(Fk[II(i), JJ(j)] ))/LL
            end
        end
    end

    M += eye(Complex128, mm)

    return factorize(M)

end



    """
    solve_smw( M, C, Linv::Function, dd1, dd2, Pm, Pp, K)
 Solves the matrix equation SMW system and computes the solution, on rectangular domains.\n
 With SMW-system matrix M, right hand side C, and matrix equation solver Linv which was used to compute M.
"""
function solve_smw( M, C::Array{Complex128,2}, Linv::Function, dd1, dd2, Pm::Function, Pp::Function, K)

    const mm::Integer = size(M,1)
    const N::Integer = sqrt(mm+4)-2      #OBS: N^2 + 4N = length(M)

    const nz::Integer = size(C,1)
    const nx::Integer = size(C,2)
    const n::Integer = nz
    const L::Integer = nz/N
    const LL::Integer = L*L

    # block extract index
    II = function (i::Integer) # z-direction, always equal
        (i-1)*L+1:i*L
    end
    JJ = function(j::Integer) #x-direction, different if boundary or interior
        ((j-3)*L+1:(j-2)*L) + 2
    end
    JJ_2 = function(j::Integer)
        (j==1)*1 + (j==2)*2 + (j==N+3)*(n+3) + (j==N+4)*(n+4)
    end

    # convert a single index k to two indeces
    # reading the matrix X by rows; left to right and top to down
    k2ij = function(k::Integer)
        j::Integer = rem(k,N+4)+(rem(k,N+4)==0)*(N+4)
        i::Integer = (k-j)/(N+4)+1
        return(i,j)
    end

    # compute the right hand side
    LinvC = Linv(C)
    b=zeros(Complex128,mm)
    for k=1:mm
        i, j = k2ij(k)
        # evaluate the linear functional
        if((j==1)||(j==2)||(j==N+3)||(j==N+4))
            b[k] = sum(sum( LinvC[II(i), JJ_2(j)] ))/L
        else
            b[k] = sum(sum( LinvC[II(i), JJ(j)] ))/LL
        end
    end
    LinvC = 0 #Clean up memory, let GC work if needed

    # solve for the coefficients
    const alpha = M\b;

    # build the solution
    Y::Array{Complex128,2} = zeros(Complex128, nz, nx);
    ek::Array{Complex128,1} = zeros(Complex128, nz);
    for k=1:mm

        i, j = k2ij(k);

        # Ek tilde implicitly created
        if (j==1)
            Y[II(i), 1] += alpha[k]*K[II(i), 1]
            ek = 0*ek
            ek[II(i)] = dd1
            Y[:, 1] += alpha[k]*Pm(ek)
        elseif (j==2)
            Y[II(i), 2] += Y[II(i),2] + alpha[k]*K[II(i), 2]
            ek = 0*ek
            ek[II(i)] = dd2
            Y[:, 1] += alpha[k]*Pm(ek)
        elseif (j==N+4)
            Y[II(i), nx] += alpha[k]*K[II(i), nx]
            ek = 0*ek
            ek[II(i)] = dd1
            Y[:, nx] += alpha[k]*Pp(ek)
        elseif (j==N+3)
            Y[II(i), nx-1] += alpha[k]*K[II(i), nx-1]
            ek = 0*ek
            ek[II(i)] = dd2
            Y[:, nx] += alpha[k]*Pp(ek)
        else
            Y[II(i),JJ(j)] += alpha[k]*K[II(i), JJ(j)]
        end

    end

    X=Linv(C-Y)
    return X

end




"""
      A wrapper that applies the preconditioner to a vector when\n
      called with the function 'solve()', which is used in GMRES
"""
    type WEP_precond_matvec
        M
        λ

        WEP_precond_matvec(nep, λ) = new(nep, λ)
    end

    function solve(A::WEP_precond_matvec, b::AbstractVector)
        return xyz
    end




