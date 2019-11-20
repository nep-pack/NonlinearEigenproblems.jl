# Types & Data structures

Nonlinear eigenvalue problems in NEP-PACK
are represented by objects of the type [`NEP`](@ref).
Each `NEP`-object needs to provide compute functions
as we describe in
[the manual page on compute functions](compute_functions.md).
Efficient compute functions are already implemented
for many common and several general types.

In the section [specific types](#Specific-types-1) below,
we list a number of common classes. As a user, first see
if your problem fits to one of those classes, as NEP-PACK has
very efficient compute functions for these classes.
If your NEP does not fit into any of the specific types, we recommend that
a user tries to specify the problem
as an [`SPMF_NEP`](@ref), which is described
in the section [general types](types.md#General-types-1).
If your problem can be phrased as a sum of two specific
or general types, it is recommended that you use the
[`SumNEP`](@ref)-type. NEP-PACK also supports efficient computation with
low-rank NEPs via the [`LowRankFactorizedNEP`](@ref).


If your NEP is not easily expressed as
an [`SPMF_NEP`](@ref), you may want to use the
[helper types](types.md#Helper-types-1).

The types also have a number of associated
operations and transformation functions.
The following example illustrates how
you can resample a NEP (by interpolation with
a Chebyshev polynomial basis in Chebyshev points
provided by the [`ChebPEP`](@ref) constructor)
and apply a NEP-solver which requires many function
evaluations, in this case [`contour_beyn`](@ref).
The two-stage solution approach is much more efficient.

```julia-repl
julia> nep_bem=nep_gallery("bem_fichera");
julia> cheb_nep=ChebPEP(nep_bem,20,0,10); # resample the NEP with 20 cheb points
 32.651313 seconds (263.16 M allocations: 36.279 GiB, 17.19% gc time)
julia> @time (λ1,v1)=contour_beyn(nep_bem,radius=[5 0.2],σ=5.0, N=100,k=10,);
180.329069 seconds (1.39 G allocations: 183.462 GiB, 13.01% gc time)
julia> @time (λ2,v2)=contour_beyn(cheb_nep,radius=[5 0.2],σ=5.0, N=100,k=10,);
  4.319376 seconds (362.34 k allocations: 8.856 GiB, 12.42% gc time)
```
Note that running the contour integral method on the
`cheb_nep` is much faster, even if we take into account that
the resampling takes some computational effort.
The computed solutions are very similar
```julia-repl
julia> λ1
2-element Array{Complex{Float64},1}:
 6.450968052414575 - 4.819767260258272e-5im
 8.105873440358572 - 0.00012794471501522612im
julia> λ2
2-element Array{Complex{Float64},1}:
 6.450968052984224 - 4.819762104884087e-5im
 8.105873439472735 - 0.0001279450670266529im
```
Moreover, if we want a very accurate solution, we can run a locally
convergence iterative method on the original problem.
It converges in very few iterations:
```julia-repl
julia> (λ2_1,v1_1)=quasinewton(nep_bem,λ=λ2[1], v=v2[:,1],logger=1);
Precomputing linsolver
iter 1 err:3.638530108313503e-12 λ=6.450968052984224 - 4.819762104884087e-5im
iter 2 err:1.2789912958165988e-14 λ=6.450968052419756 - 4.819768321350077e-5im
julia> (λ2_2,v1_2)=quasinewton(nep_bem,λ=λ2[2], v=v2[:,2],logger=1)
Precomputing linsolver
iter 1 err:3.4824421200567996e-12 λ=8.105873439472735 - 0.0001279450670266529im
iter 2 err:2.05407750614131e-14 λ=8.105873440343123 - 0.00012794469925178411im
```

!!! tip
    The use of of Chebyshev interpolation in combination with the boundary element method (but with a companion linearization approach) was presented in  [Effenberger and Kressner. "Chebyshev interpolation for nonlinear eigenvalue problems." BIT Numerical Mathematics 52.4 (2012): 933-951](https://doi.org/10.1007/s10543-012-0381-5). See also [the tutorial on boundary element method](bemtutorial.md).


## Specific types


```@docs
DEP
```

There are two types to represent PEPs natively in
NEP-PACK. You can use a monomial basis with
`PEP` or a Chebyshev basis with `ChebPEP`.

```@docs
PEP
```
```@docs
ChebPEP
```

### REP
The Rational Eigenvalue Problem is described by:

```@docs
REP
```


## General types
The basic class is the abstract class `NEP` which represents
a NEP. All other defined NEPs should inherit from `NEP`, or from a more
specialized version; see, e.g., [`ProjectableNEP`](transformations.md#NonlinearEigenproblems.NEPTypes.ProjectableNEP) or [`AbstractSPMF`](types.md#NonlinearEigenproblems.NEPTypes.AbstractSPMF).

```@docs
NEP
```


Below we list the most common types built-in to NEP-PACK, and further down how you can [access the NEP](types.md#accessNEP).
However, the structure is made for extendability, and hence it is possible for you to extend with your own class of NEPs.

## SPMF

One of the most common problem types is the `SPMF_NEP`.
SPMF is short for Sum of Products of Matrices and Functions and the NEP is described by
```math
M(λ) = \sum_{i} A_i f_i(λ).
```
The constructor of the `SPMF_NEP`, takes the
the matrices and the functions, but also a number of other (optional) parameters
which may increase performance or preserve underlying types.


```@docs
SPMF_NEP
```

```@docs
AbstractSPMF
get_Av
get_fv
```

### Projectable NEP types

There are also types associated with projection described on  [the projection manual page](innersolvers.md):
* [`ProjectableNEP`](@ref)
* [`Proj_NEP`](@ref)

### SumNEP
It is also possible to consider NEPs that are sums of other NEPs.
For such situations there are SumNEPs. Specifically `GenericSumNEP` and `SPMFSumNEP`. Both are constructed using
the function `SumNEP`.

```@docs
SumNEP
```
```@docs
GenericSumNEP
```
```@docs
SPMFSumNEP
```

### Low-rank NEPs

```@docs
LowRankFactorizedNEP
```
## CORK data types

```@docs
CORKPencil
```

```@docs
buildPencil
```
```@docs
CORKPencilLR
```

```@docs
lowRankCompress
```




## Helper types
There are also the helper types [`Mder_NEP`](@ref) and
[`Mder_Mlincomb_NEP`](@ref). These are further described in
the section about [Compute functions](compute_functions.md)
