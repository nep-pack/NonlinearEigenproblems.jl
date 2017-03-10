

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

    return K
end 

###########################################################
# Waveguide eigenvalue problem (WEP)
# Sum of products of matrices and functions (SPMF)
        
"""
 Waveguide eigenvalue problem (WEP)
Sum of products of matrices and functions (SPMF)
"""
function gallery_waveguide_spmf( )
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
    K = 1
    hx = 1
    hz = 1
    return K, hx, hz
end


    """
 Genearate the wavenumber in FD discretization for the waveguide
described by JARLEBRING.
"""
function generate_wavenumber_fd_jarlebring( nx::Integer, nz::Integer, delta)
    K = 1
    hx = 1
    hz = 1
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






