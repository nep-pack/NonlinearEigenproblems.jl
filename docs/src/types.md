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
compute_Mder(nep::NEPCore.NEP, λ::Number, i::Integer) 
compute_Mlincomb(nep::NEPCore.NEP, λ::Number, V; a)
compute_MM
```



# NEPTypes

In order to use the methods,
the user has the possibility to implement their own
problem-specific functions above, or use one of the predefined
types. The most common one is the `SPMF_NEP`.


```@docs
SPMF_NEP
```

In order to construct an `SPMF_NEP`, we need to provide
the matrices and the functions.

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






