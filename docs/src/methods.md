# NEP Methods

The NEP solver methods implemented in NEP-pack, are accessed by
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
NEPSolver.newton
NEPSolver.augnewton
NEPSolver.resinv
NEPSolver.quasinewton
NEPSolver.mslp
NEPSolver.rfi
NEPSolver.newtonqr
NEPSolver.implicitdet
```
## Projection methods
```@docs
NEPSolver.nlar
```
## Arnoldi type methods
```@docs
NEPSolver.iar
NEPSolver.tiar
NEPSolver.infbilanczos
```

## Class specific methods

### Companion linearizations

### Rational ?
