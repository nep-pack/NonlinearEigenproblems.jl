module GalleryWaveguide

using NonlinearEigenproblems.NEPCore
using NonlinearEigenproblems.NEPTypes
using NonlinearEigenproblems.LinSolvers
using NonlinearEigenproblems.Gallery
using SparseArrays
using LinearAlgebra
using Statistics
using FFTW
using IterativeSolvers
using LinearMaps

export wep_generate_preconditioner
export wep_gmres_linsolvercreator
export wep_backslash_linsolvercreator
export wep_factorized_linsolvercreator

# Specializalized NEPs
export WEP
export WEP_FD

# We overload these
import NonlinearEigenproblems.Gallery.nep_gallery
import NonlinearEigenproblems.NEPCore.compute_Mlincomb
import NonlinearEigenproblems.LinSolvers.lin_solve
import NonlinearEigenproblems.LinSolvers.default_linsolvercreator
import NonlinearEigenproblems.LinSolvers.DefaultLinSolver
import NonlinearEigenproblems.LinSolvers.BackslashLinSolver

import Base.size
import SparseArrays.issparse
import Base.*
import Base.eltype
import LinearAlgebra.ldiv!



abstract type WEP <: NEP end


include("waveguide/Waveguide.jl")
include("waveguide/waveguide_FD.jl")
include("waveguide/waveguide_FEM.jl")
include("waveguide/waveguide_preconditioner.jl")



"""
    nep_gallery(WEP, kwargs...)

Create the NEP associated with the **Waveguide Eigenvalue Problem (WEP)** found in both
* E. Ringh, and G. Mele, and J. Karlsson, and E. Jarlebring, Sylvester-based preconditioning for the waveguide eigenvalue problem, Linear Algebra Appl., 2018
* E. Jarlebring, and G. Mele, and O. Runborg, The waveguide eigenvalue problem and the tensor infinite Arnoldi method, SIAM J. Sci. Comput., 2017

# keyword parameters:
* nx::Integer = 3 * 5 * 7,         disctretization points in x-direction
* nz::Integer = 3 * 5 * 7,         disctretization points in z-direction
* benchmark_problem = 'TAUSCH',    which waveguide (TAUSCH, JARLEBRING)
* discretization::String = 'FD',   which discretization (FD)
* neptype::String = 'WEP',         NEP-format (SPMF, SPMF_PRE, WEP) later format recommended
* delta = 0.1,                     Slack from the absorbing boundary conditions
"""
function nep_gallery(::Type{T}; nx::Integer = 3*5*7, nz::Integer = 3*5*7, benchmark_problem::String = "TAUSCH", discretization::String = "FD", neptype::String = "WEP",  delta::Number = 0.1) where T<:WEP

    waveguide = uppercase(benchmark_problem)
    neptype = uppercase(neptype)
    discretization = uppercase(discretization)
    if !isodd(nz)
        error("Variable nz must be odd! You have used nz = ", nz, ".")
    end


    # Generate the matrices and formulate the problem is the sought format
    if ((neptype == "SPMF") || (neptype == "SPMF_PRE")) && (discretization == "FD")
        K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        if(neptype == "SPMF_PRE")
            pre_Schur_fact = true
        else
            pre_Schur_fact = false
        end
        nep = assemble_waveguide_spmf_fd(nx, nz, hx, Dxx, Dzz, Dz, C1, C2T, K, Km, Kp, pre_Schur_fact)

    elseif (neptype == "WEP") && (discretization == "FD")
        K, hx, hz, Km, Kp = generate_wavenumber_fd( nx, nz, waveguide, delta)
        Dxx, Dzz, Dz = generate_fd_interior_mat( nx, nz, hx, hz)
        C1, C2T = generate_fd_boundary_mat( nx, nz, hx, hz)
        nep = WEP_FD(nx, nz, hx, hz, Dxx, Dzz, Dz, C1, C2T, K, Km, Kp)

    else
        error("The NEP-type '", neptype, "' is not supported for the discretization '", discretization, "'.")
    end

    return nep
end


end
