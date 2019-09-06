# Linear solvers

Most NEP-solvers require

* [the solution of linear system of equations](linsolvers.md#Linear-system-of-equations-1), or
* [the solution of a standard eigenvalue problem](linsolvers.md#Standard-eigenvalue-problems-1).

The user can specify which linear solver or eigenvalue solver
he/she wants to use. It is also possible to use external solvers.


## Linear system of equations

Most NEP-algorithms need to solve the linear system associated with `M(λ)`.
We provide an interface to specify which solver to use or define your own solver.

```@docs
LinSolver
```
```@docs
lin_solve
```
```@docs
create_linsolver
```
```@docs
DefaultLinSolverCreator
```
```@docs
FactorizeLinSolver
```
```@docs
FactorizeLinSolverCreator
```
```@docs
BackslashLinSolver
```
```@docs
BackslashLinSolverCreator
```
```@docs
GMRESLinSolver
```
```@docs
GMRESLinSolverCreator
```


## Standard eigenvalue problems

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
