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


# EigSolvers

Some NEP-algorithms need to solve an associated linear eigenvalue problem. associated with `M(λ)`.
You will likely only need the native eigensolvers in Julia.
Nevertheless, we provide an interface to specify which solver to use or define your own solver.

```@docs
EigSolver
```
```@docs
eig_solve
```
```@docs
DefaultEigSolver
```
```@docs
NativeEigSolver
```
```@docs
NativeEigSSolver
```
