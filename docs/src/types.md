# The basic type

The basic class is the abstract class `NEP` which represents
a NEP.

```@docs
NEP
```

## Accessing the NEP

The nonlinear eigenvalue problem is defined by the data
stored in the NEP-class, and the NEP-solvers access
the data mainly through three main functions, `compute_Mder`
`compute_Mlincomb` and `compute_MM`.

```@docs
compute_Mder
```
```@docs
compute_Mlincomb
```

```@docs
compute_MM
```


# NEPTypes

In order to use the methods,
the user has the possibility to implement their own
problem-specific functions above, or use one of the predefined
types. The most common one is the `SPMF_NEP`.

## SPMF
SPMF is short for Sum of Products of Matrices and Functions and the NEP is described by
```math
M(λ) = \sum_{i} A_i f_i(λ).
```
In order to construct an `SPMF_NEP`, we need to provide
the matrices and the functions.

```@docs
SPMF_NEP
```

Many problems can be described in the class of SPMF.
There might be more specialized and efficient implementations, but on an abstract
level it may still be important to recognize the similarities.
Hence there is an abstract class `AbstractSPMF`.
```@docs
AbstractSPMF
```
All NEPs inheriting from `AbstractSPMF` needs to provide an implementation of the
following functions.
```@docs
get_Av
get_fv
```


## PEP
The Polynomial Eigenvalue Problem is described by
```math
M(λ) = \sum_{i} A_i λ^i.
```
In order to construct a `PEP`, we only need to provide
the matrices.

```@docs
PEP
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
The Rational Eigenvalue Problem is described by

```@docs
REP
```


# SumNEP
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
