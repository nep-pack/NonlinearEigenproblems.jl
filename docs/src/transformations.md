# Transforming NEPs

There are various ways to transform NEPs into other NEPs.
The simplest example is the function `shift_and_scale()`.


```@docs
shift_and_scale
```

Similarly `mobius_transform()` is more general
than `shift_and_scale` which transform
the problem using a Möbius transformation. The function `taylor_exp`
create new PEP by doing truncating a Taylor expansion.

# Projection

Several methods for NEPs are based on forming
a smaller NEP, which we will refer to as a projection:
```math
N(λ)=W^HM(λ)V,
```
where $V,W\in\mathbb{C}^{n\times p}$
and the corresponding projected problem
```math
N(λ)u=0.
```

## Types
NEPs for which this projection can be computed
inherit from `ProjectableNEP`.

```@docs
ProjectableNEP
```

The result of the
projection is represented in a `Proj_NEP`.

```@docs
Proj_NEP
```

One explicit instance is the `Proj_SPMF_NEP`.

```@docs
Proj_SPMF_NEP
```


## Associated functions

You can create a projected NEP with `create_proj_NEP`:

```@docs
create_proj_NEP
```


```@docs
set_projectmatrices!
```

```@docs
expand_projectmatrices!
```



# Deflation

Due to structure of the representation of NEPs in NEP-PACK
it is possible to do deflation, by transformation of the NEP-object.
The deflation is based on theory provided in Effenbergers thesis
and the main function consists of `deflate_eigpair`.
See also [the tutorial on deflation](deflate_tutorial.md).

```@docs
deflate_eigpair
```
