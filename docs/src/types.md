# NEPTypes

## The basic type
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

### Abstract SPMFs
Many problems can be described in the class of SPMF.
There might be more specialized and efficient implementations such as, e.g. `PEP`, `DEP` or `REP`.
However, on an abstract level it may still be important to recognize the similarities.
Hence there is an abstract class `AbstractSPMF`, which in itself inherits from [`ProjectableNEP`](transformations.md#NonlinearEigenproblems.NEPTypes.ProjectableNEP).
```@docs
AbstractSPMF
get_Av
get_fv
```


## PEP
A PEP is a matrix polynomial in `λ`:
```math
M(λ) = \sum_{i} A_i λ^i.
```
There are two types to represent PEPs natively in
NEP-PACK. You can use a monomial basis with
`PEP` or a Chebyshev basis with `ChebPEP`.

```@docs
PEP(AA::Array)
```
```@docs
ChebPEP(AA::Array)
```


## DEP
The Delay Eigenvalue Problem is described by
```math
M(λ) = -λI + \sum_{i} A_i e^{-τ_i λ}.
```

```@docs
DEP
```

## REP
The Rational Eigenvalue Problem is described by:

```@docs
REP
```
The constructor is called as:

```@docs
REP(AA,poles::Array{<:Number,1})
```


## SumNEP
It is also possible to consider NEPs that are summs of other NEPs. For such situations
there are SumNEPs. Specifically `GenericSumNEP` and `SPMFSumNEP`. Both are constructed using
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


# Accessing the NEP

The nonlinear eigenvalue problem is defined by the data
stored in the NEP-class, and the NEP-solvers access
the data mainly through three main functions, `compute_Mder`
`compute_Mlincomb` and `compute_MM`.

```@docs
compute_Mder
```

```@docs
compute_Mlincomb!
```

```@docs
compute_MM
```
