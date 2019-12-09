# Deflation

Several NEP-algorithms are able to find one eigenvalue,
but may have difficulties finding several eigenvalues.
Deflation is a transformation technique which can
transform a NEP by effectively removing computed eigenvalues,
and allowing several eigenvalues to be computed by repeated
application of the same NEP-algorithm.

NEP-PACK provides a solver-independent implementation of deflation
which can be combined (essentially) with any NEP-solver.
 NEP-PACK also has some NEP-solver deflation techniques and reconvergence avoidance techniques  incoprorated directly in the solver, e.g., in the nonlinear Arnoldi method ([`nlar`](@ref)),
the Jacobi-Davidson method ([`jd_betcke`](@ref))
and Broyden's method ([`broyden`](@ref)).

The technique takes a NEP and a solution and creates a bigger NEP
with one dimension larger but where the eigenvalue is removed from the solution set. Due to the abstraction of NEP-objects
in NEP-PACK, the deflated NEP is again a NEP and we can apply the
NEP-solver to the deflated NEP.

* Given a NEP (which can be a deflated NEP) `nep` and an eigenpair `(λ,v)` you can compute a deflated NEP by calling `dnep=`[`deflate_eigpair`](@ref)`(nep,λ,v)` and `dnep` will essentially have the same eigenvalues as `nep`, except `λ`.
* The transformation changes the eigenvectors such that the eigenvectors of `nep` and `dnep` will be different. To extract the eigenvectors (and the eigenvalues) you can call [`get_deflated_eigpairs`](@ref)`(dnep)`.

!!! note
    More elaborate deflation examples can be found in [the tutorial on deflation](deflate_tutorial.md).








## Theory

The theory follows the presentation of the technique
in [the PhD thesis of Cedric Effenberger](http://sma.epfl.ch/~anchpcommon/students/effenberger.pdf). It can be summarized as follows, in a somewhat simplified form (for the index one case).
The deflation is based on a theory for NEPs essentially stating that
if ``(s,x)`` is an eigenpair, then under certain general conditions (which we implicitly assume are satisfied), 
the extended nonlinear eigenvalue problem
```math
T(λ):=\begin{bmatrix}M(λ)&M(λ)x(s-λ)^{-1}\\ x^T & 0\end{bmatrix}
```
has the same eigenvalues as the original problem except for the eigenvalue ``s`` which is no longer part of the solution set. We have effectively removed (i.e. deflated) the eigenpair `(s,x)`. More
eigenpairs can be deflated with techniques of partial Schur
factorizations, which the user does not need to be aware of, due to
the abstraction provided by the functions below. When we create
a deflated NEP, we create the NEP ``T``.

There are several ways to represent the ``T``, which is why deflation has several
modes. If you run
```julia
julia> dnep=deflate_eigpair(nep,λ1,v1,mode=:SPMF)
```
the `dnep` will be of the type [`AbstractSPMF`](@ref). More precisely, if
```math
M(λ)=A_1f_1(λ)+\cdots+A_mf_m(λ)
```
the deflated NEP will be
```math
T(λ)=
\begin{bmatrix}A_1&0\\0 & 0\end{bmatrix}f_1(λ)+\cdots+
\begin{bmatrix}A_m&0\\0 & 0\end{bmatrix}f_m(λ)+
```
```math
\begin{bmatrix}0&A_1x\\0 & 0\end{bmatrix}\frac{f_1(λ)}{s-λ}+\cdots+
\begin{bmatrix}0&A_mx\\0 & 0\end{bmatrix}\frac{f_m(λ)}{s-λ}+
\begin{bmatrix}0&0\\x^T & 0\end{bmatrix}
```
Clearly, the deflated NEP has more SPMF-terms than the original `NEP`.
When the parameter `mode=:SPMF` is set, the deflation method will explicitly construct an [`SPMF_NEP`](@ref).
This is not recommended if you have many SPMF-terms in the original problem, but can be efficient when you
only have a few terms.
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
