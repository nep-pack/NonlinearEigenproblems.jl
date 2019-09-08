# Types & Data structures


*** TODO: Extend to describe that we have specific types, general types and composite types ***


*** Maybe add an example with `ChebPEP` to illustrate how to move between types ***


Many problems can be described in the class of SPMF.
There might be more specialized and efficient implementations such as, e.g. `PEP`, `DEP` or `REP`.
However, on an abstract level it may still be important to recognize the similarities.
Hence there is an abstract class `AbstractSPMF`, which in itself inherits from [`ProjectableNEP`](transformations.md#NonlinearEigenproblems.NEPTypes.ProjectableNEP).


The Delay Eigenvalue Problem is described by
```math
M(λ) = -λI + \sum_{i} A_i e^{-τ_i λ}.
```

A PEP is a matrix polynomial in `λ`:
```math
M(λ) = \sum_{i} A_i λ^i.
```


## Specific types


```@docs
DEP
```

There are two types to represent PEPs natively in
NEP-PACK. You can use a monomial basis with
`PEP` or a Chebyshev basis with `ChebPEP`.

```@docs
PEP(AA::Array)
```
```@docs
ChebPEP
```

### REP
The Rational Eigenvalue Problem is described by:

```@docs
REP
```
The constructor is called as:

```@docs
REP(AA,poles::Array{<:Number,1})
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
SPMF_NEP(AA::Vector{<:AbstractMatrix}, fii::Vector{<:Function};
                  Schur_fact = false,
                  check_consistency=false,
                  align_sparsity_patterns=false)
```

```@docs
AbstractSPMF
get_Av
get_fv
```

### Projectable NEP types

There are also types associated with projection described on  [the projection manual page](projection.md):
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

## Helper types
There are also the helper types [`Mder_NEP`](@ref) and
[`Mder_Mlincomb_NEP`](@ref). These are further described in
the section about [Compute functions](compute_functions.md)
