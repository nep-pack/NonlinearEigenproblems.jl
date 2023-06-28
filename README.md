# NEP-PACK

[![CI](https://github.com/nep-pack/NonlinearEigenproblems.jl/workflows/CI/badge.svg)](https://github.com/nep-pack/NonlinearEigenproblems.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/nep-pack/NonlinearEigenproblems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nep-pack/NonlinearEigenproblems.jl)


A nonlinear eigenvalue problem is the problem to determine a scalar *λ* and a vector *v* such that
*<p align="center">M(λ)v=0</p>*
where *M* is an *nxn*-matrix depending on a parameter. This package aims to provide state-of-the-art algorithms to solve this problem, as well as a framework to formulate applications and easy access to benchmark problems. This currently includes (but is not restricted to) Newton-type methods, Subspace methods, Krylov methods, contour integral methods, block methods, companion matrix approaches. Problem transformation techniques such as scaling, shifting, deflating are also natively supported by the package.


# How to use it?

On Julia 1.X, install it as a registered package by typing `] add ...` at the REPL-prompt:
```
julia> ]
(v1.0) pkg> add NonlinearEigenproblems
```

After that, check out "Getting started" in

<p align="center"><a href="https://nep-pack.github.io/NonlinearEigenproblems.jl">NEP-PACK online user's guide</a></p>

or read the preprint: https://arxiv.org/abs/1811.09592

## GIT Version

If you want the cutting edge development version and not the latest release, install it with the URL:
```
julia> ]
(v1.0) pkg> add git://github.com/nep-pack/NonlinearEigenproblems.jl.git
```
## NEP solvers

Features and solvers (see documentation https://nep-pack.github.io/NonlinearEigenproblems.jl/methods/ for further information and references):

* Arnoldi/Krylov type
    * NLEIGS
    * Infinite Arnoldi method: (iar)
    * Tensor infinite Arnoldi method  (tiar)
    * Infinite bi-Lanczos (infbilanczos)
    * Infinite Lanczos (ilan)
    * AAA CORK (AAAeigs)
* Projection methods
    * Jacobi-Davidson (jd_effenberger)
    * Jacobi-Davidson (jd_betcke)
    * Nonlinear Arnoldi method (nlar)
    * Common Rayleigh-Ritz projection interface
* Contour integral methods
    * Beyn's contour integral method
    * Block SS (Higher moments) contour integral method of Asakura & Sakurai
    * Common quadrature interface for parallelization
* Newton & Rayleigh type:
    * Classical Newton-Raphson
    * Augmented Newton
    * Residual inverse iteration
    * Quasi-Newton
    * Block Newton
    * Rayleigh functional iteration (RFI a, b)
    * Newton-QR
    * Implicit determinant method
    * Broyden's method
* Companion matrices
    * First companion form
    * Companion form for Chebyshev polynomials
* Other
    * Chebyshev interpolation
    * Transformations
    * Rayleigh-Ritz (`ProjNEP` and `inner_solve`)
    * Problem gallery (including access to the NLEVP-gallery)
    * Deflation (Effenberger style)


# Development

Core developers (alphabetical): Max Bennedich, Elias Jarlebring (www.math.kth.se/~eliasj), Giampaolo Mele (www.math.kth.se/~gmele), Emil Ringh (www.math.kth.se/~eringh), Parikshit Upadhyaya (https://www.kth.se/profile/pup/). Thanks to A Koskela for involvement in initial version of the software.

# How to cite

If you find this software useful please cite

```bibtex
@Misc{,
  author = 	 {E. Jarlebring and M. Bennedich and G. Mele and E. Ringh and P. Upadhyaya},
  title = 	 {{NEP-PACK}: A {Julia} package for nonlinear eigenproblems},
  year = 	 {2018},
  note = 	 {https://github.com/nep-pack},
  eprint = {arXiv:1811.09592},
}
```
If you use a specific method, please also give credit to the algorithm researcher.
Reference to a corresponding algorithm paper can be found by in, e.g., by writing `?resinv`.

Some links below are developer info on KTH. We will migrate them soon:


* NEP-page style "guide": https://github.com/nep-pack/NonlinearEigenproblems.jl/wiki/Style-guidelines-and-notes

* GIT-workflow: https://github.com/nep-pack/NonlinearEigenproblems.jl/wiki/Git-workflow

* GIT-usage @ KTH: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki

* NEP-methods @ KTH: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/NEP-methods

* NEP-applications @ KTH: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/Applications
