# Transformations.

There are various ways to transform NEPs into other NEPs.
The simplest example is the function `shift_and_scale()`.


```@docs
shift_and_scale
```

Similarly `mobius_transform()` is more general
than `shift_and_scale` which transform
the problem using a Möbius transformation. The function `taylor_exp`
create new PEP by doing truncating a Taylor expansion.


```@docs
mobios_transform
```

** TODO: This should be renamed, e.g. taylor_exp **

```@docs
transform_to_pep
```
