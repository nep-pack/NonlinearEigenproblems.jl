# NEP-PACK

[![Build Status](https://img.shields.io/travis/nep-pack/NonlinearEigenproblems.jl.svg)](https://travis-ci.org/nep-pack/NonlinearEigenproblems.jl)
[![Coverage Status](https://img.shields.io/coveralls/github/nep-pack/NonlinearEigenproblems.jl.svg?label=coveralls)](https://coveralls.io/github/nep-pack/NonlinearEigenproblems.jl?branch=master)
[![codecov](https://img.shields.io/codecov/c/github/nep-pack/NonlinearEigenproblems.jl.svg?label=codecov)](https://codecov.io/gh/nep-pack/NonlinearEigenproblems.jl)

A nonlinear eigenvalue problem is the problem to determine a scalar *λ* and a vector *v* such that
*<p align="center">M(λ)v=0</p>*
where *M* is an *nxn*-matrix depending on a parameter. This package aims to provide state-of-the-art algorithms to solve this problem, as well as a framework to formulate applications and easy access to benchmark problems. This currently includes (but are not restricted to) Newton-type methods, Subspace methods, Krylov methods, contour integral methods, block methods, companion. Problem transformation techniques such as scaling, shifting, deflating are also directly supported.  

# How to use it?


We are currently in a transition phase to Julia 1.0 (or julia 0.7). Stay tuned.

On julia 0.6.4 you can install it by running
```
julia> Pkg.clone("git://github.com/nep-pack/NonlinearEigenproblems.jl.git");
```
Due to the transition, you need also switch to git-tag v0.1. 

After that, check out "Getting started" in the user documentation
https://nep-pack.github.io/NonlinearEigenproblems.jl


# Development

The main work of NEP-PACK has been done in a closed repository at KTH, but as of May 2018 the development is carried out in a public github repo.

Core developers (alphabetical): Max Bennedich, Elias Jarlebring (www.math.kth.se/~eliasj), Giampaolo Mele (www.math.kth.se/~gmele), Emil Ringh (www.math.kth.se/~eringh), Parikshit Upadhyaya (https://www.kth.se/profile/pup/). Thanks to A Koskela for involvement in initial version of the software.

Developers who want to push directly to the repository can use the following installation
```
julia> # Pkg.rm("NonlinearEigenproblems")  #
julia> Pkg.clone("git@github.com:nep-pack/NonlinearEigenproblems.jl.git");
```

# How to cite

If you find this software useful please cite

```bibtex
@Misc{,
  author = 	 {E. Jarlebring, G. Mele, E. Ringh, P. Upadhyaya},
  title = 	 {NEP-PACK: A nonlinear eigenvalue problem in Julia},
  year = 	 {2018},
  note = 	 {https://github.com/nep-pack},
}
```
If you use a specific method, please also give credit to the algorithm researcher.
Reference to a corresponding algorithm paper can be found by in, e.g., by writing `?resinv`.

Some links below are developer info on KTH. We will migrate them soon:


* Checklist for first public version: https://github.com/nep-pack/NonlinearEigenproblems.jl/issues/26

* NEP-page style "guide": https://github.com/nep-pack/NonlinearEigenproblems.jl/wiki/Style-guidelines-and-notes

* GIT-usage: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki

* NEP-methods: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/NEP-methods

* NEP-applications: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/Applications
