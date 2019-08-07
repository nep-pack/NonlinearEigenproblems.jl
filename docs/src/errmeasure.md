# Measuring the error

All iterative algorithms need some form of termination
criteria. In NEP-PACK, all NEP-solvers provide
the possibility to specify the desired tolerance,
as well as how the error is measured or estimated.
The tolerance is specified in the kwarg  `tol` (which is a real number)
and the way to measure the error is given in `errmeasure`.

## Standard usage

NEP-PACK comes with several ways to measure errors for many NEP-types.

* `errmeasure=ResidualErrmeasure`: The error is estimated by the use of the residual norm:
```math
\mathrm{err}=\frac{\|M(λ)v\|}{\|v\|}.
```
* `errmeasure=BackwardErrmeasure`: The error is estimated by using the backward error bounds. This error measure will not work for all NEPs. An implementation is provided for any `AbstractSPMF`. If your NEP is an `AbstractSPMF` with terms:
  ```math
  M(λ)=A_1f_1(λ)+\cdots+A_mf_m(λ)
  ```
  the error will be estimated by
  ```math
  \mathrm{err}=\frac{\|M(λ)v\|}{\|v\|}\frac{1}{\|A_1\|_F|f_1(λ)|+\cdots+\|A_m\|_F|f_m(λ)|}.
  ```
  In other words, the `BackwardErrmeasure` is a weighting of the `ResidualErrmeasure`.
* `errmeasure=DefaultErrmeasure`: When this `errmeasure` is specified, NEP-PACK tries to determine a error measure for you. In general, `BackwardErrmeasure` will be preferred if possible. This behavior may change in future versions of NEP-PACK.

Example: Most NEP-solvers take the `errmeasure` as an kwarg.
```julia-repl
julia> nep=nep_gallery("qdep0");
julia> # Solve the problem to residual norm 1e-8
julia> (λ,v)=mslp(nep,errmeasure=ResidualErrmeasure,tol=1e-8)
julia> norm(compute_Mlincomb(nep,λ,v))/norm(v) # It's smaller than tol?
3.503700819937386e-9
julia> nep isa AbstractSPMF # Is it an AbstractSPMF so we can use BackwardErrmeasure?
true
julia> (λ,v)=mslp(nep,errmeasure=BackwardErrmeasure,tol=1e-10)
julia> factor=abs(fv[1](λ))*norm(Av[1])+
     abs(fv[2](λ))*norm(Av[2])+abs(fv[3](λ))*norm(Av[3]);
julia> norm(compute_Mlincomb(nep,λ,v))/(norm(v)*factor)
1.659169482386331e-11
```

## User defined error measure

There are two ways that a user can specify how to measure the error.

### User defined error 1: Function handle

The user can provide a function handle
which is called to obtain the error. The `errmeasure` can be a function,
which takes two parameters as input `(λ,v)` and returns
the error (or estimate thereof).

The most common situation is that you want to report the
error (as a function of iteration) with a reference solutions.
If we want to get
a very accurate approximation of the true error, we can run the
algorithm twice, and the second time we run the algorithm
we use the result of the first run as a reference.

```julia-repl
julia> nep=nep_gallery("qdep0");
julia> v0=ones(size(nep,1));
julia> (λref,v)=resinv(nep,v=v0,λ=230^2+1im,logger=1);
julia> myerrmeasure = (λ,v) -> abs(λ-λref)/abs(λ);
julia> (λ,v)=resinv(nep,v=v0,λ=250^2+1im,logger=1,tol=1e-10,errmeasure=myerrmeasure);
Iteration:  1 errmeasure:1.274091618457501296e-01
Iteration:  2 errmeasure:9.535794095609478882e-01
...
Iteration: 49 errmeasure:1.269396691930517923e-10
Iteration: 50 errmeasure:7.608430406801784718e-11
```

### User defined error 2: A user defined type

Due to the multiple dispatch and handling of types in Julia, code
may run faster if one uses types instead of function handles. It is
possible to do the same simulation as above with a user defined
type.

You first need to define a new type
```julia-repl
julia> struct RefErrmeasure <: Errmeasure; nep::NEP; end
```
The error measure should then provided in the function
`estimate_error`:
```julia-repl
julia> v0=ones(size(nep,1));
julia> (λref,v)=resinv(nep,v=v0,λ=230^2+1im,logger=1);
julia> function NonlinearEigenproblems.estimate_error(e::RefErrmeasure,λ,v)
         return abs(λ-λref)/abs(λ);
       end
julia> (λ,v)=resinv(nep,v=v0,λ=250^2+1im,logger=1,tol=1e-10,errmeasure=RefErrmeasure);
Iteration:  1 errmeasure:1.274091618457501296e-01
...
Iteration: 49 errmeasure:1.269396691930517923e-10
Iteration: 50 errmeasure:7.608430406801784718e-11
```


## As a NEP-solver developer

NEP-solvers should use the `Errmeasure` as follows. The NEP-solver should take
as input a `ErrmeasureType`. This corresponds to either a function or a type, but if you follow the procedure below, you will not have to worry about that.

```julia
function mysolver(nep::NEP;errmeasure::ErrmeasureType=DefaultErrmeasure)
```

Before the main iteration, you need to initialize the error measure
computation. The precomptued data
is stored in a variable typically called `ermdata`:
```julia
   ermdata=init_errmeasure(errmeasure,nep);
```

In the main for loop you want to call the `estimate_error` function:

```julia
for k=1:maxit
    err=estimate_error(ermdata,λ,v)
    if (err < 1e-10)
       return (λ,v)
    end
    ....

end
```

## Methods and types

```@docs
Errmeasure
```
```@docs
estimate_error
```
```@docs
init_errmeasure
```
```@docs
 DefaultErrmeasure
```
```@docs
ResidualErrmeasure
```
```@docs
BackwardErrmeasure
```
