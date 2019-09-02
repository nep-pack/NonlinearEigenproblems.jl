# Inner solvers

Many NEP-solvers are based on a computation of a
 projected problem, i.e., if ``V,W\in\mathbb{R}^{n\times p}``
we need to solve the (smaller) NEP
```math
W^HM(λ)Vz=0
```
This is sometimes called a nonlinear Rayleigh-Ritz procedure,
or a direct projection.

NEP-PACK provides a framework to handle projected problems.
You can in principle use any of the NEP-solvers to
solve a projected problem. As a user, this
is specified in the `inner_solver_method`
keyword argument.

If you wish to use the infinite Arnoldi method
to handle the project solves in the nonlinear
Arnoldi method, you can do the following:

```julia-repl
julia> nep=nep_gallery("dep0_tridiag");
julia> λ,v=nlar(nep,neigs=1,inner_solver_method=IARInnerSolver(),logger=1);
Using inner solver IARInnerSolver(1.0e-13, 80, :ones, false, NonlinearEigenproblems.NEPSolver.iar)
iter 1 err:0.05095382004494062 λ=0.7579134426195271 - 0.03707164055891316im
iter 2 err:0.00031997693290503965 λ=-0.00010049358638757657 + 0.0001763732030940319im
iter 3 err:6.563177508431498e-6 λ=-0.0005335154073888051 - 4.498082881742902e-6im
iter 4 err:8.037612383023366e-9 λ=-0.0005259618179586685 + 2.806438064753968e-9im
iter 5 err:3.386041718599221e-11 λ=-0.0005259683064505526 + 1.3119844096548209e-11im
iter 6 err:3.4499779767886924e-13 λ=-0.0005259682927702825 - 1.0404412824030558e-13im
iter 7 err:5.365662696809372e-15 λ=-0.0005259682929434785 - 2.0247938561528697e-17im
****** 1 converged to eigenvalue: -0.0005259682929434785 - 2.0247938561528697e-17im errmeasure:5.365662696809372e-15
```
The logging of the inner solver is controlled by the kwarg `inner_logger`.
This produces very verbose output illustrating
also the convergence of the inner solve:
```julia-repl
julia> λ,v=nlar(nep,neigs=1,inner_solver_method=IARInnerSolver(),logger=1,inner_logger=1);
Using inner solver IARInnerSolver(1.0e-13, 80, :ones, false, NonlinearEigenproblems.NEPSolver.iar)
-
--
---
----
-----
------
-------
--------
---------
----------
-----------
------------
=------------
+-------------
iter 1 err:0.04860520921206162 λ=0.8437634284420165 + 0.005742468974957178im
-
--
---
+---
iter 2 err:0.0008771218464072076 λ=-0.00065275394814732 - 0.0008601482370586537im
-
--
---
...
```
## Inner solvers

The inner solvers inherit from `InnerSolver`.
The following inner solvers are available by default.


```@docs
NewtonInnerSolver
```

```@docs
IARInnerSolver
```

```@docs
IARChebshevInnerSolver
```

```@docs
ContourBeynInnerSolver
```

```@docs
PolyeigInnerSolver
```

```@docs
SGIterInnerSolver
```


## Advanced usage

You can define your own inner solver by
inheriting from `InnerSolver` and implementing
the function `inner_solve`. Since the `inner_solve`
obtains information from the solver via
keyword arguments, you need to end your
method signature with `kwargs...)`.

```@docs
InnerSolver
```

```@docs
inner_solve
```

## Rayleigh functional computation

```@docs
compute_rf
```