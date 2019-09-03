# Deflation

Due to structure of the representation of NEPs in NEP-PACK
it is possible to do deflation, by transformation of the NEP-object.
The deflation is based on theory provided in Effenbergers thesis
and the main function consists of [`deflate_eigpair`](@ref).
See also [the tutorial on deflation](deflate_tutorial.md).

** TODO: This description needs to be extended (maybe move theory from tutorial) **


```@docs
deflate_eigpair
```

```@docs
get_deflated_eigpairs
```