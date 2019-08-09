# NEP Methods

The NEP solver methods implemented in NEP-PACK, are accessed by
the functions below. The functions all return ``λ,v,w`` where
``λ`` is either a number (eigenvalue) a vector of eigenvalues
``v`` is either a vector containing an eigenvector
or a matrix whose columns corresponding to the eigenvectors.

The first parameter optional parameter in all NEP solver methods
is a type. This type specifies which arithmetic should be used
for the algorithm.

Example:

```julia-repl
julia> nep=nep_gallery("dep0")
julia> λ,v=augnewton(Complex128,nep,v=ones(5))
(0.8347353572199425 + 0.0im, Complex{Float64}[0.480386+0.0im, 0.0631636+0.0im, -0.136405+0.0im, 0.214274+0.0im, 0.378581+0.0im])
julia> typeof(λ)
Complex{Float64}
julia> λ,v=augnewton(Float16,nep,v=ones(5))
(Float16(0.8223), Float16[0.47388, 0.063904, -0.13843, 0.21692, 0.38306])
julia> typeof(λ)
Float16
```


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


## Class specific methods

### Companion linearizations
```@docs
companion
polyeig
```

### Rational ?
