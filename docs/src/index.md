
# NEPCore Documentation



```@docs
NEP
```

```@docs
SPMF_NEP
```

```@docs
SPMF_NEP(AA::Array, fii::Array{Function,1}) 
```

```@docs
compute_Mder
```


```@docs
compute_Mlincomb
```

Let's try an equation $x=f(x)$. \(a \ne 0\)

```math
x=x_1+1
```

# Compiling the documentation

Compile this documentation site by running:
```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ julia --color=yes make.jl
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ mkdocs build --clean
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ firefox site/index.html
```

More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)




# NEP methods

## Newton type methods
```@docs
NEPSolver.newton
NEPSolver.resinv
NEPSolver.quasinewton
NEPSolver.mslp
```
## Projection methods


