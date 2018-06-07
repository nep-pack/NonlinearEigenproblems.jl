# NEP-PACK

A nonlinear eigenvalue problem is the problem to determine a scalar λ and a vector v such that
<p align="center">M(λ)v=0</p>
where M is an nxn-matrix depending on a parameter. This package is a first version of methods for nonlinear eigenvalue problems in julia. 

# How to use it?

It's easy. Installation in a one-liner:
```
julia> Pkg.clone("git://github.com/nep-pack/NonlinearEigenproblems.jl.git");
```
Check out "Getting started" in the user documentation
https://nep-pack.github.io/NonlinearEigenproblems.jl

# Developers

Core developers (alphabetical): Elias Jarlebring (www.math.kth.se/~eliasj), Giampaolo Mele (www.math.kth.se/~gmele), Emil Ringh (www.math.kth.se/~eringh), Parikshit Upadhyaya (https://www.kth.se/profile/pup/). Thanks to A Koskela for involvement in initial version of the software.

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
If you use a specific method, please also give credit to the method developer.
Reference to the method developer can be found by in, e.g., by writing `?resinv`.

Links below are developer info on KTH. We will migrate them soon: 

Documentation page: https://gitr.sys.kth.se/pages/nep-pack/nep-pack-alpha/


* Checklist for first public version: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/issues/26

* NEP-page style "guide": https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/Style-guidelines-and-notes

* GIT-usage: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki

* NEP-methods: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/NEP-methods

* NEP-applications: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/Applications


