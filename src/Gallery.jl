  """
  Module containing a gallery of examples of nonlinear eigenvalue problems.\\
  Look at the function 'nep_gallery()' for further instructions.
  """
module Gallery
    using NonlinearEigenproblems.Serialization
    using ..NEPCore
    using ..NEPTypes
    using Random
    using LinearAlgebra
    using SparseArrays
    using PolynomialRoots
    using SpecialFunctions

    export nep_gallery

    function __init__()
        # Add the search-path to the extra galleries - so that unloaded modules can be easily loaded from top level
        this_path = string(@__DIR__, "/gallery_extra")
        if !(this_path in LOAD_PATH)
            push!(LOAD_PATH, this_path)
        end
    end

    # Add the search-path to the extra galleries - temporarily so so needed files can be included
    push!(LOAD_PATH, string(@__DIR__, "/gallery_extra"))
    include("gallery_extra/basic_random_examples.jl")
    include("gallery_extra/gallery_examples.jl")
    include("gallery_extra/distributed_example.jl")
    include("gallery_extra/periodic_dde.jl")
    include("gallery_extra/NLEVP_native.jl")

  """
     nep=nep_gallery(name)
     nep=nep_gallery(name,params)
     nep=nep_gallery(name,params;kwargs)

**Collection of nonlinear eigenvalue problems.**\\
Returns a NEP object from a gallery of examples of nonlinear eigenvalue problems.\\
The parameter `name` decides which NEP.

# Supported problems:
The following list describes the NEP with a certain `name` and the associated parameters (`params`) and keyword arguments (`kwargs`), if any.\\
\\
\\
     `dep0`\\
Create a random delay eiganvalue problem with one delay tau = 1
* one optional parameter determining the size (default = 5)
\\
     `dep0_sparse`\\
Create a random delay eiganvalue problem with sparse matrices and one delay tau = 1
* two optional parameter determining the size (default = 5) and the fill (default = 0.25)
\\
      `dep0_tridiag`\\
Create a random delay eiganvalue problem with sparse tridiaognal matrices and one delay tau = 1
* one optional parameter determining the size (default = 100)
\\
      `dep_symm_double`\\
Create delay eiganvalue problem with double eigenvalues and sparse symmetric matrices and one delay tau = 1\\
Examle from H. Voss and M. M. Betcke, Restarting iterative projection methods for Hermitian nonlinear eigenvalue problems with minmax property, Numer. Math., 2017
* one optional parameter determining the size (default = 100)
\\
     `dep_double`\\
Create problem with a double non-semisimple eigenvalue in λ=3πi\\
      Example from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012\\
\\
     `dep1`\\
A delay eigenvalue problem with one eigenvalue equal to one.\\
\\
     `pep0`\\
Create a random polynomial eigenvalue problem
* one optional parameter determining the size (default = 200)
 \\
      `pep0_sym`\\
Creates a random symmetric polynomial eigenvalue problem
* one optional parameter determining the size (default = 200)
\\
     `pep0_sparse`\\
Creates a random polynomial eigenvalue problem with sparse matrices
* two optional parameters determining the size (default = 200) and the fill (default = 0.03)
\\
     `real_quadratic`\\
Creates a quadratic problem with real eigenvalues\\
          Four smallest eigenvalues of the problem:\\
          -2051.741417993845\\
          -182.101627437811\\
          -39.344930222838\\
          -4.039879577113\\
\\
     `dep_distributed`\\
Creates the NEP associated with example in E. Jarlebring and W. Michiels and K. Meerbergen,  The infinite Arnoldi method and an application to time-delay systems with distributed delays,\\
Delay Systems - Methods, Applications and New Trends, 2012\\
         Some correct eigenvalues:\\
         -0.400236388049641 + 0.970633098237807i\\
         -0.400236388049641 - 0.970633098237807i\\
          2.726146249832675 + 0.000000000000000i\\
         -1.955643591177653 + 3.364550574688863i\\
         -1.955643591177653 - 3.364550574688863i\\
          4.493937056300693 + 0.000000000000000i\\
         -1.631513006819252 + 4.555484848248613i\\
         -1.631513006819252 - 4.555484848248613i\\
         -1.677320660400946 + 7.496870451838560i\\
         -1.677320660400946 - 7.496870451838560i\\
\\
     `qdep0` \\
Quadratic delay eigenvalue problem in S. W. Gaaf and E. Jarlebring, The infinite Bi-Lanczos method for nonlinear eigenvalue problems, SIAM J. Sci. Comput., 2017 \\
\\

     `qdep1` \\
Quadratic delay eigenvalue problem in E. Jarlebring and W. Michiels and K. Meerbergen,  A linear eigenvalue algorithm for the  nonlinear eigenvalue problem, Numer. Math., 2011 \\
\\
     `qep_fixed_eig`\\
A quadratic eigenvalue problem with chosen eigenvalues
* two optional parameters determining the size (default = 5)
       and a vector containing the eigenvalues (default = randn)
\\
     `neuron0`\\
A DEP that stems from L. P. Shayer and S. A. Campbell, Stability, bifurcation and multistability in a system of two coupled neurons with multiple time delays,
    SIAM J. Applied Mathematics, 2000. It is also a benchmark example in DDE-BIFTOOL\\
\\
     `beam`\\
The DEP modelling a beam with delayed stabilizing feedback described in R. Van Beeumen, E. Jarlebring, and W. Michiels, A rank-exploiting infinite Arnoldi algorithm for nonlinear eigenvalue problems, 2016.\\
The A1-term has rank one.
* one optional parameter which is the size of the matrix (defalut = 100)
\\
     `sine` \\
The NEP formed by the sum of a polynomial and a sine-function in "A rank-exploiting infinite Arnoldi algorithm for nonlinear eigenvalue problems", R. Van Beeumen, E. Jarlebring and W. Michiels, 2016. The sine-term has rank one.\\
\\
     `nlevp_native_gun`\\
The benchmark problem from the NLEVP-collection called "gun", represented in the native NEP-PACK format. B.-S. Liao, Z. Bai, L.-Q. Lee, and K. Ko. Nonlinear Rayleigh-Ritz iterative method for solving large scale nonlinear eigenvalue pro blems.  Taiwan. Journal of Mathematics, 14(3):869–883, 2010\\
\\
     `nlevp_native_fiber`\\
The benchmark problem from the NLEVP-collection called "fiber", represented in the native NEP-PACK format. One of terms in this problem is approximated by interpolation, and may not always coincide with the benchmark. Kaufman, L. 2006. Eigenvalue problems in fiber optic design. SIAM J. Matrix Anal. Appl. 28, 1, 105–117.  and Huang, X., Bai, Z., and Su, Y. 2010. Nonlinear rank-one modification of the symmetric eigenvalue problem. J. Comput. Math. 28, 2, 218–234.\\

# Example
```julia-repl
julia> nep=nep_gallery("dep0",100);
julia> norm(compute_Mlincomb(nep,1.0+1.0im,ones(size(nep,1))))
104.76153002802755
```

# See also the following galleries:
* GalleryNLEVP
* GalleryWaveguide
  """
    nep_gallery(name::String, params...; kwargs...)=nep_gallery(NEP, name, params...; kwargs...)
    function nep_gallery(::Type{T}, name::String, params...; kwargs...) where T<:NEP
        global Gallery
        if !haskey(GALLERY, name)
            error("$name not supported")
        end

        return GALLERY[name](params...; kwargs...)
    end

    const GALLERY = Dict(
        "dep0" => dep0,
        "dep1" => dep1,
        "dep0_sparse" => dep0_sparse,
        "dep0_tridiag" => dep0_tridiag,
        "dep_symm_double" => dep_symm_double,
        "dep_double" => dep_double,
        "pep0" => pep0,
        "pep0_sym" => pep0_sym,
        "pep0_sparse" => pep0_sparse,
        "real_quadratic" => real_quadratic,
        "dep_distributed" => gallery_dep_distributed,
        "qdep0" => qdep0,
        "qdep1" => qdep1,
        "qep_fixed_eig" => qep_fixed_eig,
        "periodicdde" => (params...; kwargs...) -> periodic_dde_gallery(PeriodicDDE_NEP, params...; kwargs...),
        "neuron0" => neuron0,
        "nlevp_native_gun" => nlevp_native_gun,
        "nlevp_native_fiber" => nlevp_native_fiber,
        "beam" => beam,
        "sine" => sine_nep,
    )


end # end of Gallery-module
