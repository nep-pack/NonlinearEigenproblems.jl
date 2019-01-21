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
    include("gallery_extra/bem_hardcoded/bem_hardcoded.jl");
    include("gallery_extra/load_dtn_dimer.jl");

  """
     nep=nep_gallery(name)
     nep=nep_gallery(name,params)
     nep=nep_gallery(name,params;kwargs)

Collection of nonlinear eigenvalue problems.
Returns a NEP object from a gallery of examples of nonlinear eigenvalue problems.
The parameter `name` decides which NEP.

# Supported problems:
The following list describes the NEP with a certain `name` and the associated parameters (`params`) and keyword arguments (`kwargs`), if any.\\

* `dep0`
    Create a random delay eiganvalue problem with one delay tau = 1.\\
    One optional `params` determining the size (default = 5)
* `dep0_sparse`\\
    Create a random delay eiganvalue problem with sparse matrices and one delay tau = 1.\\
    Two optional `params` determining the size (default = 5) and the fill (default = 0.25)

* `dep0_tridiag`\\
    Create a random delay eiganvalue problem with sparse tridiaognal matrices and one delay tau = 1.\\
    One optional `params` determining the size (default = 100)

* `dep_symm_double`\\
    Create delay eiganvalue problem with double eigenvalues and sparse symmetric matrices and one delay tau = 1.\\
    Examle from H. Voss and M. M. Betcke, Restarting iterative projection methods for Hermitian nonlinear eigenvalue problems with minmax property, Numer. Math., 2017\\
    One optional `params` determining the size (default = 100)

* `dep_double`\\
    Create problem with a double non-semisimple eigenvalue in λ=3πi.\\
    Example from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012

* `dep1`\\
    A delay eigenvalue problem with one eigenvalue equal to one.

* `pep0`\\
    Create a random polynomial eigenvalue problem.\\
    One optional `params` determining the size (default = 200)

* `pep0_sym`\\
    Creates a random symmetric polynomial eigenvalue problem.\\
    One optional `params` determining the size (default = 200)

* `pep0_sparse`\\
    Creates a random polynomial eigenvalue problem with sparse matrices.\\
    Two optional `params` determining the size (default = 200) and the fill (default = 0.03)

* `real_quadratic`\\
    Creates a quadratic problem with real eigenvalues.\\
          Four smallest eigenvalues of the problem:\\
          -2051.741417993845\\
          -182.101627437811\\
          -39.344930222838\\
          -4.039879577113

* `dep_distributed`\\
    Creates the NEP associated with example in E. Jarlebring and W. Michiels and K. Meerbergen,
    The infinite Arnoldi method and an application to time-delay systems with distributed delays,
    Delay Systems - Methods, Applications and New Trends, 2012.\\
         Some correct eigenvalues:
         `-0.400236388049641 + 0.970633098237807i`,\\
         `-0.400236388049641 - 0.970633098237807i`,\\
         ` 2.726146249832675 + 0.000000000000000i`,\\
         `-1.955643591177653 + 3.364550574688863i`,\\
         `-1.955643591177653 - 3.364550574688863i`,\\
         ` 4.493937056300693 + 0.000000000000000i`,\\
         `-1.631513006819252 + 4.555484848248613i`,\\
         `-1.631513006819252 - 4.555484848248613i`,\\
         `-1.677320660400946 + 7.496870451838560i`,\\
         `-1.677320660400946 - 7.496870451838560i`

* `qdep0` \\
    Quadratic delay eigenvalue problem in S. W. Gaaf and E. Jarlebring, The infinite Bi-Lanczos method for
    nonlinear eigenvalue problems, SIAM J. Sci. Comput., 2017

* `qdep1` \\
    Quadratic delay eigenvalue problem in E. Jarlebring and W. Michiels and K. Meerbergen,
    A linear eigenvalue algorithm for the  nonlinear eigenvalue problem, Numer. Math., 2011

* `qep_fixed_eig`\\
    A quadratic eigenvalue problem with chosen eigenvalues.\\
    Two optional `params` determining the size (default = 5)
    and a vector containing the eigenvalues (default = randn)

* `neuron0`\\
    A DEP that stems from L. P. Shayer and S. A. Campbell, Stability, bifurcation and multistability
    in a system of two coupled neurons with multiple time delays,
    SIAM J. Applied Mathematics, 2000. It is also a benchmark example in DDE-BIFTOOL

* `schrodinger_movebc` \\
    This NEP stems from the discretization of a Schrödinger equation as described in the NEP-PACK online tutorial. The nonlinearity contains ``sinh()``, ``cosh()`` and ``sqrt()``. The optional parameters are size of discretization `n`  and domain and potential description `L0`,`L1`,`α` and `V0`.

* `beam`\\
    The DEP modelling a beam with delayed stabilizing feedback described in R. Van Beeumen, E. Jarlebring, and W. Michiels,
    A rank-exploiting infinite Arnoldi algorithm for nonlinear eigenvalue problems, 2016.\\
    The A1-term has rank one.\\
    One optional `params` which is the size of the matrix (defalut = 100)

* `sine` \\
    The NEP formed by the sum of a polynomial and a sine-function in "A rank-exploiting infinite Arnoldi
    algorithm for nonlinear eigenvalue problems", R. Van Beeumen, E. Jarlebring and W. Michiels, 2016. The sine-term has rank one.

    The MATLAB-package "NLEVP: A Collection of Nonlinear Eigenvalue Problems, ACM Transactions on Mathematical Software 39(2), January 2011,
    T. Betcke, N. J. Higham, V. Mehrmann, Ch. Schröder, F. Tisseur" provides a number of benchmark problems for NEPs.
    These are available in NEP-PACK in two different ways. We have native implementations of some problems (referred to as `nlevp_native_`)
    and the separate `GalleryNLEVP`. The native implementation is preferred since the `GalleryNLEVP`
    interfaces with MATLAB and is therefore considerably slower.

* `nlevp_native_gun`\\
    The benchmark problem from the NLEVP-collection called "gun", represented in the native NEP-PACK format.
    B.-S. Liao, Z. Bai, L.-Q. Lee, and K. Ko. Nonlinear Rayleigh-Ritz iterative method for solving large scale
    nonlinear eigenvalue problems.  Taiwan. Journal of Mathematics, 14(3):869–883, 2010

* `nlevp_native_cd_player`\\
    The benchmark problem from the NLEVP-collection called "cd_player", represented in the native NEP-PACK format.
    Y. Chahlaoui, and P. M. Van Dooren, Benchmark examples for model reduction of linear time-
    invariant dynamical systems. In Dimension Reduction of Large-Scale Systems, P. Benner, V. Mehrmann,
    and D. C. Sorensen, Eds. Lecture Notes in Computational Science and Engineering Series, vol. 45.
    Springer-Verlag, Berlin, 380–392, 2005.\\
    and\\
    P. M. R. Wortelboer, M. Steinbuch, and  O. H. Bosgra, Closed-loop balanced reduction with
    application to a compact disc mechanism. In Selected Topics in Identification, Modeling and Control.
    Vol. 9. Delft University Press, 47–58, 1996.\\
    and\\
    W. Draijer, M. Steinbuch, and  O. H. Bosgra, Adaptive control of the radial servo system of a
    compact disc player. Automatica 28, 3, 455–462. 1992.\\


* `nlevp_native_fiber`\\
    The benchmark problem from the NLEVP-collection called "fiber", represented in the native NEP-PACK format.
    One of terms in this problem is approximated by interpolation, and may not always coincide with the benchmark.
    L. Kaufman, Eigenvalue problems in fiber optic design. SIAM J. Matrix Anal. Appl. 28, 1, 105–117, 2006.\\
    and\\
    X. Huang, Z. Bai, and Y. Su, Nonlinear rank-one modification of the symmetric eigenvalue problem. J. Comput. Math. 28, 2, 218–234, 2010


* `nlevp_native_hadeler`\\
    The benchmark problem from the NLEVP-collection called "hadeler", represented in the native NEP-PACK format. The problem is of the form ``M(λ)=(e^λ-1)B+A0+A2λ^2``. \\
    Hadeler K.  P.  1967.  Mehrparametrige  und  nichtlineare  Eigenwertaufgaben. Arch.  Rational  Mech. Anal. 27, 4, 306–328.\\

* `nlevp_native_pdde_stability`\\
    The benchmark problem from the NLEVP-collection called "pdde_stability", represented in the native NEP-PACK format.
    This problem is a quadratic eigenvalue with arbitrary given size `n`. See
    E. Jarlebring, The Spectrum of Delay-Differential Equations:
    Numerical Methods, Stability and Perturbation, PhD thesis,
    TU Braunschweig, Institut Computational Mathematics, Germany, 2008 and \\
    H. Fassbender, N. Mackey, D. S. Mackey and C. Schroeder, Structured
    Polynomial Eigenproblems Related to Time-Delay Systems, ETNA, 2008, vol 31, pp 306-330\\

* `nlevp_native_loaded_string`\\
    The benchmark problem from the NLEVP-collection called "pdde_stability", represented in the native NEP-PACK format. The parameters are (n,kappa,m) where n is the size, and the NEP is a SPMF with rational terms and the coefficient matrices are rank one modifications of Toeplitz matrices.\\
   S. I. Solov"ev. Preconditioned iterative methods for a class of nonlinear eigenvalue problems. Linear Algebra Appl., 415 (2006), pp.210-229.

* `bem_fichera`\\
    Represents a boundary element discretization of Helmholtz equation for a domain consisting of the unit cube, except one removed corner (Fichera corner). The mesh is hardcoded. The parameter N determines the size of the problem (default N = 5). The model stems from the model in these papers:\\
   Steinlechner, A boundary element method for solving PDE eigenvalue problems, M. Steinlechner, bachelor thesis, ETH Zürich, 2010\\
   Effenberger and Kressner, Chebyshev interpolation for nonlinear eigenvalue problems, BIT Numerical Mathematics, December 2012, Volume 52, Issue 4, pp 933–951

* `dtn_dimer`\\
   NEP described in J. Araujo-Cabarcas, C. Engström and E. Jarlebring, Efficient resonance computations for Helmholtz problems based on a Dirichlet-to-Neumann map, J. Comput. Appl. Math., 330:177-192, 2018  (http://arxiv.org/pdf/1606.09547) where the nonlinearity stems from the application of Dirichlet-to-Neumann map. The problem contains quotients of Bessel functions and derivatives of Bessel functions. This NEP takes two parameters: `data_dir::String` and `l::Int`. The `data_dir` specifies the directory of the dowloaded FEM-matrices (available here https://umu.app.box.com/s/b52yux3z9rcl8y0l7la22k0vi062cvu5). The integer l specifies the number of DtN-terms: 2*l+1.   \\
   J. Araujo-Cabarcas, C. Engström and E. Jarlebring, Efficient resonance computations for Helmholtz problems based on a Dirichlet-to-Neumann map, J. Comput. Appl. Math., 330:177-192, 2018  (http://arxiv.org/pdf/1606.09547).


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
        "schrodinger_movebc" => schrodinger_movebc,
        "neuron0" => neuron0,
        "nlevp_native_gun" => nlevp_native_gun,
        "nlevp_native_hadeler" => nlevp_native_hadeler,
        "nlevp_native_cd_player" => nlevp_native_cd_player,
        "nlevp_native_fiber" => nlevp_native_fiber,
        "nlevp_native_pdde_stability" => nlevp_native_pdde_stability,
        "nlevp_native_loaded_string" => nlevp_native_loaded_string,
        "beam" => beam,
        "sine" => sine_nep,
        "bem_fichera" => bem_fichera,
        "dtn_dimer" => load_dtn_dimer
    )


end # end of Gallery-module
