# Deflation

Due to structure of the representation of NEPs in NEP-PACK
it is possible to do deflation, by transformation of the NEP-object.
The deflation is based on theory provided in Effenbergers thesis
and the main function consists of [`deflate_eigpair`](@ref).
See also [the tutorial on deflation](deflate_tutorial.md).

** TODO: This description needs to be extended (maybe move theory from tutorial) **


## The theory in the background

The deflation is based on a theory for NEP essentially stating that
if ``(s,x)`` is an eigenpair, then the extended nonlinear eigenvalue problem
```math
T(λ):=\begin{bmatrix}M(λ)&M(λ)x(s-λ)^{-1}\\ x^T & 0\end{bmatrix}
```
has the same eigenvalues as the original problem (under certain quite general
conditions which are assumed to be satisfied). More
eigenpairs can be deflated with techniques of partial Schur
factorizations (which the user does not need to use). When we create
a deflated NEP, we create the NEP `T`.

There are several ways to represent the ``T``, which is why deflation has several
modes. If you run
```julia
julia> dnep=deflate_eigpair(nep,λ1,v1,mode=:SPMF)
```
the `dnep` will be of the type `AbstractSPMF`. More precisely, if
```math
M(λ)=A_1f_1(λ)+\cdots+A_mf_m(λ)
```
the deflated NEP will be
```math
T(λ)=
\begin{bmatrix}A_1&0\\0 & 0\end{bmatrix}f_1(λ)+\cdots+
\begin{bmatrix}A_m&0\\0 & 0\end{bmatrix}f_m(λ)+
\begin{bmatrix}0&A_1x\\0 & 0\end{bmatrix}\frac{f_1(λ)}{s-λ}+\cdots+
\begin{bmatrix}0&A_mx\\0 & 0\end{bmatrix}\frac{f_m(λ)}{s-λ}+
\begin{bmatrix}0&0\\x^T & 0\end{bmatrix}
```
Clearly, the deflated NEP will have more SPMF-terms, and
the `mode=:SPMF`, is not recommended if you have many SPMF-terms.
(Some additional exploitation is however implemented, since we can use
the fact that the introduced terms are of low rank, and
therefore naturally represented as a `LowRankFactorizedNEP`.)

If you select `mode=:Generic`, the compute functions are implemented
without the use of SPMF, and can be more efficient
if the NEP has many SPMF-terms.
When `mode=:MM` the compute-functions are all implemented
by calls to `compute_MM`. This will not be efficient if
`compute_Mder(nep,λ,der)` where  `der>0` is needed.

## Functions

```@docs
deflate_eigpair
```

```@docs
get_deflated_eigpairs
```



![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_DEFLATIONM)
