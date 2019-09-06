# Projection

Many NEP-solvers are based on a computation of a solution
to a projected problem, i.e., if ``V,W\in\mathbb{R}^{n\times p}``
we need to solve the (smaller) NEP
```math
W^HM(λ)Vz=0
```
This is sometimes called a nonlinear Rayleigh-Ritz procedure,
or a direct projection. These are *inner solvers* for many NEP-solvers.


NEP-PACK provides a framework to handle projected problems
and inner solves. This is implemented
into two separate components:

* [Projection](@ref): As a user (or NEP-solver developer) you can create a new object corresponding to the projection. In NEP-PACK, the projection is again an object of with type inheriting from [`NEP`](@ref). More precisely, it is a [`Proj_NEP`](@ref) which you normally create with the function [`create_proj_NEP`](@ref).
* [Inner solvers](@ref): Since the projected problem is again a `NEP`, in principle any of the NEP-solvers of this package can be used. This is handled by the `InnerSolver` objects which are wrappers for corresponding NEP-solvers such that we can pass appropriate parameters to the inner soler. The inner solver is controlled by the `inner_solver_method` keyword in many NEP-solvers. By default [`DefaultInnerSolver`](@ref) is used.



As a NEP-user, you often do not need to care about how the
projection is handled, e.g., if you use the type [`SPMF_NEP`](@ref)
with only a few terms. For instance,
if you wish to use the infinite Arnoldi method ([`iar`](@ref))
to handle the project solves in the nonlinear
Arnoldi method ([`nlar`](@ref)), you can do the following:

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

The logging of the inner solver is controlled by the kwarg `inner_logger`,
which follows the same framework as the standard [NEP-PACK Logger](logger.md).
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

[Rayleigh functional computation](@ref), which corresponds to projection
with $p=1$, is also handled with this framework.

## Inner solvers

The inner solvers inherit from [`InnerSolver`](@ref).
The following inner solvers are available by default.


```@docs
NewtonInnerSolver
```

```@docs
IARInnerSolver
```

```@docs
IARChebInnerSolver
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

```@docs
NleigsInnerSolver
```

```@docs
DefaultInnerSolver
```


## Inner solvers: Advanced usage

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


## Projection

The NEP-PACK functionality for projected problems
are represented by [projection types](@ref).
Normally, the projection is created by
[`create_proj_NEP`](@ref) from a standard NEP.
After creating a projected NEP, you can set
the projection subspace (represented by the
matrices `V` and `W`) using
[`set_projectmatrices!`](@ref) or
[`expand_projectmatrices!`](@ref).
```julia-repl
julia> A=[1 0 0; 0 1.0 0; 0 0 1]; B=[1 2 3; 3 3 3 ; 4 -1 -1.0];
julia> nep=SPMF_NEP([A, B], [s->s, s->s^5]);
julia> pnep=create_proj_NEP(nep);
julia> W=[4 1 ; 6 1  ; 6.0 2]; V=[3 3;3 4.0;4.0 -1];
julia> set_projectmatrices!(pnep,W,V); # modifies pnep
julia> λ=3.0+1im;
julia> W'*compute_Mder(nep,λ)*V
2×2 Array{Complex{Float64},2}:
 -3366.0+92958.0im  -2238.0+61334.0im
  -690.0+19290.0im   -513.0+13909.0im
julia> compute_Mder(pnep,λ)
2×2 Array{Complex{Float64},2}:
 -3366.0+92958.0im  -2238.0+61334.0im
  -690.0+19290.0im   -513.0+13909.0im
```
Effectively, the `Proj_NEP` creates [compute functions](compute_functions.md),
which are designed to be as efficient as possible.

### Projection functions

You can create a projected NEP with `create_proj_NEP`:

```@docs
create_proj_NEP
```


```@docs
set_projectmatrices!
```

```@docs
expand_projectmatrices!
```

### Projection types
NEPs for which this projection can be computed
inherit from `ProjectableNEP`.

```@docs
ProjectableNEP
```

The result of the
projection is represented in a `Proj_NEP`.

```@docs
Proj_NEP
```

One explicit instance is the `Proj_SPMF_NEP`.

```@docs
Proj_SPMF_NEP
```




## Rayleigh functional computation


```@docs
compute_rf
```
