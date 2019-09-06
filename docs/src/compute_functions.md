# Compute functions


The nonlinear eigenvalue problems in NEP-PACK are defined
by the data stored in the corresponding NEP-class.
The advised way NEP-solvers access the data is to do it
through three main functions,
which take the NEP-object as input.
* [`compute_Mder`](@ref): Computes a given derivative of the matrix function $M(λ)$.
* [`compute_Mlincomb!`](@ref) (or [`compute_Mlincomb`](@ref)): Computes a linear combination of derivatives $M(λ)$
* [`compute_MM`](@ref): Computes the block residual.

The choice of these functions as the fundamental way
to access a NEP is a balancing between what applications
can provide and NEP-solvers need.

A user who needs a new class of NEPs (which is
not available in among the standard types)
is advised to use the helper functions
[`Mder_NEP`](@ref) and/or
[`Mder_Mlincomb_NEP`](@ref) rather than
reimplementing the compute-functions, since
the helper types are more user friendly.
Implementation of your own NEP-type is only
advised if needed for efficiency reasons.


As a NEP-solver developer,
`compute_Mlincomb`-calls are preferred over
`compute_Mder`-calls, for the same reasons that
algorithms that only require matrix vector products can be
easier to use in a given application than an iterative
algorithm using only matrix vector products. It is in general
also more efficient although they produce the same result
up to round-off errors:
```julia-repl
julia> using BenchmarkTools;
julia> n=1000; p=10;
julia> nep=DEP([randn(n,n), randn(n,n)];
julia> V=randn(n,p);
julia> @btime compute_Mlincomb(nep,1.0,V);
  478.316 μs (19 allocations: 80.78 KiB)
julia> @btime for k=1:p; z[:]+=compute_Mder(nep,1.0,k)*V[:,k]; end
  78.510 ms (183 allocations: 465.71 MiB)
```
The `compute_Mlincomb`-function exist in two variants,
where `compute_Mlincomb!` may modify the `V`-matrix,
but in general require less memory allocations.

For a type where only `compute_Mder` is implemented,
the `compute_Mlincomb`-functionality can be provided
by delegating using the function
[`compute_Mlincomb_from_Mder`](@ref), such that
methods which require `compute_Mlincomb` can be used.


## Compute-functions documentation
```@docs
compute_Mder
```

```@docs
compute_Mlincomb
```
```@docs
compute_Mlincomb!
```

```@docs
compute_MM
```


## Type helpers


```@docs
Mder_NEP
Mder_Mlincomb_NEP
```

```@docs
compute_Mlincomb_from_Mder
```
```@docs
compute_Mlincomb_from_MM
```
