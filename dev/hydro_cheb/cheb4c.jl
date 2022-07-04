using LinearAlgebra

function cheb4c(N)

    # The function x,D4 = cheb4c(N) computes the fourth derivative
    # matrix D4 on Chebyshev interior points, incorporating the
    # clamped boundary conditions.
    # It is adapted from the function with the same name in dmsuite
    # developed for Matlab by JAC Weideman and SC Reddy.
    #
    # Adapted by M. Beneitez at KTH Mechanics.
    # beneitez@mech.kth.se
    #
    # Input:
    # N-2:  Order of the differentiation matrix.
    #
    # Output:
    # yF: Interior Chebyshev points
    # D4: Size [N-2,N-2] Fourth derivative matrix
    #

    eye = Matrix{Float64}(I, N-2, N-2);

    n1 = Int(floor(N/2-1)); n2  = Int(ceil(N/2-1));     # Indices used for flipping trick.

    k = collect(1:N-2);                        # Compute theta vector.
    th = k.*pi./(N-1)

    x = sin.(pi.*collect(N-3:-2:3-N)./(2*(N-1))); # Compute Chebyshev points.

    s = [sin.(th[1:n1]);reverse(sin.(th[1:n2]), dims = 1)]

    alpha = s.^4;
    beta1 = -4*s.^2 .*x./alpha;
    beta2 = 4*(3*x.^2 .-1)./alpha;
    beta3 = 24*x./alpha;
    beta4 = 24 ./alpha;
    B = [beta1'; beta2'; beta3'; beta4'];

    T = repeat(th/2,1,N-2);
    DX = 2*sin.(T'+T).*sin.(T'-T);          # Trigonometric identity. 
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:], dims = 2), dims=1)];   # Flipping trick. 
    DX[eye.==true] = ones(N-2,1);                      # Put 1's on the main diagonal of DX.

    ss = s.^2 .*(-1.0).^k;
    S  = repeat(ss, 1, N-2);
    C  = S./S';

    Z  = 1 ./DX;
    Z[eye.==true] = zeros(N-2,1);

    X  = Z';
    X  = X[X .!= 0.0];
    X  = reshape(X,N-3,N-2);

    Y  = ones(N-3,N-2);
    D  = Matrix{Float64}(I, N-2, N-2);
    DM = zeros(N-2,N-2,4);

    for ell = 1:4
        Y = cumsum([B[ell,:]'; ell.*Y[1:N-3,:].*X], dims = 1);
        D = ell.*Z.*(C.*repeat(diag(D), 1, N-2)-D);
        D[eye.==true] = Y[N-2,:];
        DM[:,:,ell] = D;
    end

    D4 = DM[:,:,4];

    return x,D4
    
    
end
