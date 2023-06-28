# Linear solvers

Most NEP-solvers require

* [the solution of linear system of equations](linsolvers.md#Linear-system-of-equations-1), or
* [the solution of a standard eigenvalue problem](linsolvers.md#Standard-eigenvalue-problems-1).

The user can specify which linear solver or eigenvalue solver
he/she wants to use. It is also possible to use external or user-defined
solvers.


## Linear system of equations


As a user, you can provide a creator object to many NEP-solvers
via the keyword argument `linsolvercreator`.
The creator object corresponds to one (or several) linear
system solvers. By default, [`DefaultLinSolverCreator`](@ref)
is used which tries to determine an appropriate linear
solver based on the NEP-type. In the next section,
we list the default linear solvers.


If you wish to define your own linear solver,
you need to define your own type inheriting from [`LinSolver`](@ref)
as well as a `LinSolverCreator`.
See the documenation for [`LinSolver`](@ref) and
[the tutorial on linear solvers](tutorial_linsolve.md).

### `LinSolver`-objects and creators

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
```@docs
DefaultLinSolverCreator
```
```@docs
DeflatedNEPLinSolver
```
```@docs
DeflatedNEPLinSolverCreator
```

### Advanced usage

```@docs
LinSolver
```
```@docs
lin_solve
```
```@docs
create_linsolver
```


## Standard eigenvalue problems

Some NEP-algorithms need to solve an associated linear eigenvalue problem,
associated with `M(λ)`.
We provide the possibility to
use Julia-native eigenvalue solvers, and
an interface which allows you to define your own solver.
By default, [`DefaultEigSolver`](@ref) is specified, which
tries to determine an appropriate eigenvalue solver
based on the NEP-type.

!!! tip
    The NEP-solvers [`mslp`](@ref) and [`sgiter`](@ref) require eigenvalue solvers and take the keyword argument `eigsolver`.

```@docs
DefaultEigSolver
```
```@docs
EigenEigSolver
```
```@docs
ArnoldiEigSolver
```
```@docs
EigSolver
```
```@docs
eig_solve
```
