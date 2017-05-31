
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


