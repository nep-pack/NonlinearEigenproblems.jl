

export gallery_waveguide

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
function gallery_waveguide( nx::Integer = 3*5*7, nz::Integer = 3*5*7, waveguide::String = "TAUSCH", discretization::String = "FD", NEP_format_type::String = "SPMF",  delta = 0.1)
    waveguide = uppercase(waveguide)
    NEP_format_type = uppercase(NEP_format_type)
    discretization = uppercase(discretization)


    # Generate the matrices for the sought waveguide
    if discretization == "FD"
        K, hx, hz = generate_wavenumber_fd( nx, nz, waveguide, delta)
       # Dxx, Dzz, Dz, C1, C2T = generate_fd_matrices( nx, nz, hx, hz)
    elseif discretization == "FEM"
        error("FEM discretization of WEP is not implemented yet.")
        #K = generate_wavenumber_fem( nx, nz, waveguide, delta)
    else
        error("The discretization '", discretization, "' is not supported.")
    end


    # Formulate the problem is the sought format
    if NEP_format_type == "SPMF"
        nep = K #assemble_waveguide_spmf TODO: This is a Placeholder!
    else
        error("The NEP-format '", NEP_format_type, "' is not supported.")
    end


    return nep
end 

###########################################################
# Waveguide eigenvalue problem (WEP)
# Sum of products of matrices and functions (SPMF)
        
"""
 Waveguide eigenvalue problem (WEP)
Sum of products of matrices and functions (SPMF)
"""
function assemble_waveguide_spmf( )
    #return SPMF_NEP(AA,fii::Array)
end    


###########################################################
# Generate Wavenumber FINITE DIFFERENCE
    """
 Genearate a wavenumber for Finite Difference.
"""
function generate_wavenumber_fd( nx::Integer, nz::Integer, wg::String, delta)
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
function generate_wavenumber_fd_tausch( nx::Integer, nz::Integer, delta)
    xm = 0;
    xp = 2/pi + 0.4;
    zm = 0;
    zp = 1;

    xm = xm - delta;
    xp = xp + delta;

    # Domain (First generate including the boundary)
    X = collect(linspace(xm, xp, nx+2));
    Z = collect(linspace(zm, zp, nz+1));
    # Removing the boundary
    X = X[2:end-1];
    Z = Z[2:end];

    hx = X[2]-X[1];
    hz = Z[2]-Z[1];

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

    return K, hx, hz
end


    """
 Genearate the wavenumber in FD discretization for the waveguide
described by JARLEBRING.
"""
function generate_wavenumber_fd_jarlebring( nx::Integer, nz::Integer, delta)
    xm = -1;
    xp = 1;
    zm = 0;
    zp = 1;

    xm = xm - delta;
    xp = xp + delta;

    # Domain (First generate including the boundary)
    X = collect(linspace(xm, xp, nx+2));
    Z = collect(linspace(zm, zp, nz+1));
    # Removing the boundary
    X = X[2:end-1];
    Z = Z[2:end];

    hx = X[2]-X[1];
    hz = Z[2]-Z[1];

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
            k3 *(x.>-1) .* (x.<=(-1+1)) .* (z.>0.5) .* (z.*x_ones-(x.*z_ones)/2.<=1) +
            k2 *(x.>-1) .* (x.<=(-1+1)) .* (z.>0.5) .* (z.*x_ones-(x.*z_ones)/2.>1) +
            k3 *(x.>-1) .* (x.<=(-1+1)) .* (z.<=0.5) .* (z.*x_ones+(x.*z_ones)/2.>0) +
            k2 *(x.>-1) .* (x.<=(-1+1)) .* (z.<=0.5) .* (z.*x_ones+(x.*z_ones)/2.<=0);
        end

    K = k(X', Z).^2;
    return K, hx, hz
end


###########################################################
# Generate Wavenumber FINITE ELEMENT
    """
 Genearate a wavenumber for Finite Element.
"""
function generate_wavenumber_fem( nx::Integer, nz::Integer, wg::String, delta)
    # Use early-bailout principle. If a supported waveguide is found, compute and return. Otherwise end up here and throw and error
    error("No wavenumber loaded: The given Waveguide '", wg ,"' is not supported in 'FEM' discretization.")
end






