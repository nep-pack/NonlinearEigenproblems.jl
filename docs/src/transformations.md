# Transformations

Due to the object oriented way of handle NEPs in NEP-PACK,
a NEP-object can be transformed to another NEP-object
in a number of ways. There is support
for:

* [variable transformations](transformations.md#Change-of-variables-1),
* [expansions](transformations.md#Expansions-1), and
* [deflation](transformations.md#Deflation-1).


## Change of variables

```@docs
shift_and_scale
```

```@docs
mobius_transform
```
## Expansions

```@docs
taylor_exp_pep
```

```@docs
interpolate
```

## Deflation

A NEP can be transformed to another NEP by extending
the problem in a way that it essentially removes
eigenvalues. This type of deflation is described
on [the manual page for deflation](deflation.md).