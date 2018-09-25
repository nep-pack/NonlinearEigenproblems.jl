# LinSolvers

Most NEP-algorithms need to solve the linear system associated with `M(λ)`.
We provide an interface to specify which solver to use or define your own solver.

```@docs
LinSolver
```
```@docs
lin_solve
```
```@docs
DefaultLinSolver
```
```@docs
default_linsolvercreator
```
```@docs
BackslashLinSolver
```
```@docs
backslash_linsolvercreator
```
```@docs
GMRESLinSolver
```
```@docs
gmres_linsolvercreator
```
