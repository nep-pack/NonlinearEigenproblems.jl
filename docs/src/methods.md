# NEP-Solvers

The NEP solver methods implemented in NEP-PACK, are accessed by
the functions below. The functions all return ``λ,v`` where
``λ`` is either a number (eigenvalue) a vector of eigenvalues
``v`` is either a vector containing an eigenvector
or a matrix whose columns corresponding to the eigenvectors.
Two-sided methods may return ``λ,v,w`` where ``w`` are the left eigenvectors.

The first optional parameter in all NEP solver methods
is a type. This type specifies which arithmetic should be used
for the algorithm.

Example:

```julia-repl
julia> nep=nep_gallery("dep0");
julia> λ,v=augnewton(ComplexF64,nep,v=ones(5))
(-0.15955391823299256 + 0.0im, Complex{Float64}[0.12505315954062152 + 0.0im, 0.8475907515488971 + 0.0im, -0.10910413290558324 + 0.0im, 0.027714719799174125 + 0.0im, 0.10874550201689052 + 0.0im])
julia> typeof(λ)
Complex{Float64}
julia> λ,v=augnewton(Float16,nep,v=ones(5))
(Float16(-0.718), Float16[0.435, 0.6606, -0.205, -0.1445, 0.254])
julia> typeof(λ)
Float16
```

The NEP-solvers can be separated into the following types (with some overlap):

* [Newton type methods](methods.md#Newton-type-methods-1)
* [Projection methods](methods.md#Projection-methods-1)
* [Contour integral methods](methods.md#Contour-integral-methods-1)
* [Arnoldi and Krylov based methods](methods.md#Arnoldi-and-Krylov-based-methods-1)
* [Class specific methods](methods.md#Class-specific-methods-1)

## Newton type methods
```@docs
newton
```
```@docs
augnewton
```
```@docs
resinv
```
```@docs
quasinewton
```
```@docs
mslp
```
```@docs
sgiter
```
```@docs
rfi
```
```@docs
rfi_b
```
```@docs
blocknewton
```
```@docs
newtonqr
```
```@docs
implicitdet
```
```@docs
broyden
```

## Projection methods
```@docs
nlar
jd_betcke
jd_effenberger
```
The following NEP-solvers can also be seen as
projection methods:
* [`iar`](@ref), [`tiar`](@ref), [`iar_chebyshev`](@ref),
* [`nleigs`](@ref).

## Contour integral methods

```@docs
contour_beyn
contour_block_SS
```

## Arnoldi and Krylov based methods



### IAR
The Infinite ARnoldi method.
```@docs
iar
```

### IAR Chebyshev
A Chebyshev version of the IAR method.

```@docs
iar_chebyshev
```
For the `iar_chebyshev` the following `compute_y0_cheb` method is needed, in order
to avoid explicit conversions between the Chebyshev basis and the monimial basis.

```@docs
NEPSolver.compute_y0_cheb
```


### TIAR
The Tensor Infinite ARnoldi method.

```@docs
tiar
```

### Infinite Lanczos based methods
The Infinite Bi-Lanczos method.
```@docs
infbilanczos
```
The Infinite Lanczos method, for symmetric NEPs
```@docs
ilan
```


### NLEIGS
```@docs
nleigs
```

### AAA-EIGS
```@docs
AAAeigs
```

## Class specific methods

### Companion linearizations
```@docs
companion
polyeig
```
