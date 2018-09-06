

###########################################################
# Generate discretization matrices for FINITE DIFFERENCE
    """
 Genearate the discretization matrices for the interior for Finite Difference.
"""
function generate_fd_interior_mat( nx, nz, hx, hz)
    ex = ones(nx)
    ez = ones(nz)

    # DISCRETIZATION OF THE SECOND DERIVATIVE
    Dxx = spdiagm(-1 => ex[1:end-1], 0 => -2*ex, 1 => ex[1:end-1])
    Dzz = spdiagm(-1 => ez[1:end-1], 0 => -2*ez, 1 => ez[1:end-1])
    #IMPOSE PERIODICITY IN Z-DIRECTION
    Dzz[1, end] = 1;
    Dzz[end, 1] = 1;

    Dxx = Dxx/(hx^2);
    Dzz = Dzz/(hz^2);

    # DISCRETIZATION OF THE FIRST DERIVATIVE
    Dz  = spdiagm(-1 => -ez[1:end-1], 1 => ez[1:end-1])

    #IMPOSE PERIODICITY
    Dz[1, end] = -1;
    Dz[end, 1] = 1;

    Dz = Dz/(2*hz);

    return Dxx, Dzz, Dz
end


    """
 Genearate the discretization matrices for the boundary for Finite Difference.
"""
function generate_fd_boundary_mat( nx, nz, hx, hz)
    # BUILD THE SECOND BLOCK C1
    e1 = spzeros(nx,1)
    e1[1] = 1
    en = spzeros(nx,1)
    en[end] = 1
    Iz = sparse(1.0I, nz, nz)
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
    X = range(xm, stop = xp, length = nx+2)
    hx = step(X);
    X = collect(X);
    Z = range(zm, stop = zp, length = nz+1)
    hz = step(Z);
    Z = collect(Z);
    # Removing the boundary
    X = X[2:end-1];
    Z = Z[2:end];


    # The actual wavenumber
    k1 = sqrt(2.3)*pi;
    k2 = sqrt(3)*pi;
    k3 = pi;
    k = function(x,z)
            z_ones = ones(size(z)) #Placeholder of ones to expand x-vector
            k1*(x .<= 0) .* z_ones +
            k2*(x.>0) .* (x.<=2/pi) .* z_ones +
            k2*(x.>2/pi) .* (x.<=(2/pi+0.4)) .* (z.>0.5) +
            k3*(x.>2/pi) .* (z.<=0.5) .* (x.<=(2/pi+0.4)) +
            k3*(x.>(2/pi+0.4)) .* z_ones;
        end

    K = k(X', Z).^2;
    Km = k(-Inf, 1/2)[1];
    Kp = k(Inf, 1/2)[1];
    return K, hx, hz, Km, Kp
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
    X = range(xm, stop = xp, length = nx+2)
    hx = step(X);
    X = collect(X);
    Z = range(zm, stop = zp, length = nz+1)
    hz = step(Z);
    Z = collect(Z);
    # Removing the boundary
    X = X[2:end-1];
    Z = Z[2:end];


    # The actual wavenumber
    k1 = sqrt(2.3)*pi;
    k2 = 2*sqrt(3)*pi;
    k3 = 4*sqrt(3)*pi;
    k4 = pi;
    k = function(x,z)
            z_ones = ones(size(z)) #Placeholder of ones to expand x-vector
            x_ones = ones(size(x)) #Placeholder of ones to expand z-vector
            k1 *(x.<=-1) .* z_ones  +
            k4 *(x.>1) .* z_ones  +
            k4 *(x.>(-1+1.5)) .* (x.<=1) .* (z.<=0.4) +
            k3 *(x.>(-1+1)) .* (x.<=(-1+1.5)) .* z_ones +
            k3 *(x.>(-1+1.5)) .* (x.<=1) .* (z.>0.4) +
            k3 *(x.>-1) .* (x.<=(-1+1)) .* (z.>0.5) .* (z.*x_ones-(x.*z_ones)/2 .<= 1) +
            k2 *(x.>-1) .* (x.<=(-1+1)) .* (z.>0.5) .* (z.*x_ones-(x.*z_ones)/2 .> 1) +
            k3 *(x.>-1) .* (x.<=(-1+1)) .* (z.<=0.5) .* (z.*x_ones+(x.*z_ones)/2 .> 0) +
            k2 *(x.>-1) .* (x.<=(-1+1)) .* (z.<=0.5) .* (z.*x_ones+(x.*z_ones)/2 .<= 0);
        end

    K = k(X', Z).^2;
    Km = k(-Inf, 1/2)[1];
    Kp = k(Inf, 1/2)[1];
    return K, hx, hz, Km, Kp
end
