using LinearAlgebra
using ToeplitzMatrices

function chebdif(N, M)

    # The function x,DM = chebdif(N,M) computes the differentiation
    # matrices D1, D2, ..., DM on Chebyshev nodes.
    # It is adapted from the function with the same name in dmsuite
    # developed for Matlab by JAC Weideman and SC Reddy.
    #
    # Adapted by M. Beneitez at KTH Mechanics.
    # beneitez@mech.kth.se
    #
    # Input:
    # N:  Size of the differentiation matrix.
    # M:  Number of derivates required
    # Note: 0 < M <= N-1.
    #
    # Output:
    # DM: DM[1:N,1:N,ell] contains the ell-th derivative matrix, ell=1...M
    #

    eye = Matrix{Float64}(I, N, N);        # Identity matrix.

    n1 = Int(floor(N/2)); n2  = Int(ceil(N/2));     # Indices used for flipping trick.

    k = collect(0:N-1);                        # Compute theta vector.
    th = k.*pi./(N-1)

    x = sin.(pi.*collect(N-1:-2:1-N)./(2*(N-1))); # Compute Chebyshev points.

    T = repeat(th/2,1,N);                
    DX = 2*sin.(T'+T).*sin.(T'-T);          # Trigonometric identity. 
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:], dims = 2), dims=1)];   # Flipping trick. 
    DX[eye.==true] = ones(N,1);                      # Put 1's on the main diagonal of DX.

    C = Matrix(Toeplitz((-1).^vec(k),(-1).^vec(k)));            # C is the matrix with 
    C[1,:] = C[1,:].*2; C[N,:] = C[N,:].*2;     # entries c[k]/c[j]
    C[:,1] = C[:,1]/2; C[:,N] = C[:,N]/2

    Z = 1 ./DX;                              # Z contains entries 1/(x[k]-x[j])  
    Z[eye.==true] = zeros(N,1);              # with zeros on the diagonal.

    D  = eye;                                 # D contains diff(). matrices.
    DM = zeros(N,N,M);

    for ell = 1:M
        D = ell*Z.*(C.*repeat(diag(D),1,N) - D); # Off-diagonals
        D[eye.==true] = -sum(D', dims = 1);                            # Correct main diagonal of D
        DM[:,:,ell] = D;                                   # Store current D in DM
    end

    return x,DM
    
end
