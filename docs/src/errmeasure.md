# Measuring the error

All iterative algorithms need some form of termination
criteria. In NEP-PACK, all NEP-solvers provide
the possibility to specify the desired tolerance,
as well as how the error is measured or estimated.
The tolerance is specified in the kwarg  `tol` (which is a real number)
and the way to measure the error is given in `errmeasure`.

## Standard usage

NEP-PACK comes with several ways to measure errors for many NEP-types.

* `errmeasure=`[`ResidualErrmeasure`](@ref)`(nep)`: The error is estimated by the use of the residual norm:
```math
\mathrm{err}=\frac{\|M(λ)v\|}{\|v\|}.
```
* `errmeasure=`[`StandardSPMFErrmeasure`](@ref)`(nep)`: The error is estimated by using backward error theory. This error measure will not work for all NEPs. An implementation is provided for any `AbstractSPMF`. If your NEP is an `AbstractSPMF` with terms:
  ```math
  M(λ)=A_1f_1(λ)+\cdots+A_mf_m(λ)
  ```
  the error will be estimated by
  ```math
  \mathrm{err}=\frac{\|M(λ)v\|}{\|v\|}\frac{1}{\|A_1\|_F|f_1(λ)|+\cdots+\|A_m\|_F|f_m(λ)|}.
  ```
  In other words, the `StandardSPMFErrmeasure` is a weighting of
  the `ResidualErrmeasure`.
* `errmeasure=`[`DefaultErrmeasure`](@ref)`(nep)`: When this `errmeasure` is specified, NEP-PACK tries to determine a error measure for you. In general, `StandardSPMFErrmeasure` will be preferred if possible. This behavior may change in future versions of NEP-PACK.

* `errmeasure=`[`EigvalReferenceErrmeasure`](@ref)`(nep,λref)`: This errmeasure is used when an exact (or very accurate) eigenvalue is already known. Typically, if you wish to visualize the eigenvalue error of a specific method, you run the method twice and use the result of the first run as to instantiate this error measure and get real eigenvalue errors as output.

* `errmeasure=(λ,v)-> compute_error(λ,v)`: A user defined error measure can be specified using a function. The function should be take an eigenpair as input, and return a real value. See [`ErrmeasureType`](@ref) for an example.

Example: Most NEP-solvers take the `errmeasure` as an kwarg.
```julia-repl
julia> nep=nep_gallery("qdep0");
julia> # Solve the problem to residual norm 1e-8
julia> (λ,v)=mslp(nep,errmeasure=ResidualErrmeasure(nep),tol=1e-8)
julia> norm(compute_Mlincomb(nep,λ,v))/norm(v) # It's smaller than tol?
3.503700819937386e-9
julia> nep isa AbstractSPMF # Is it an AbstractSPMF so we can use StandardSPMFErrmeasure?
true
julia> (λ,v)=mslp(nep,errmeasure=StandardSPMFErrmeasure(nep),tol=1e-10)
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
julia> (λref,_)=resinv(nep,v=v0,λ=-0.1,logger=1);
julia> myerrmeasure = (λ,v) -> abs(λ-λref)/abs(λ);
julia> (λ,v)=resinv(nep,v=v0,λ=-0.1,logger=1,tol=1e-10,errmeasure=myerrmeasure);
Precomputing linsolver
iter 1 err:0.02854168838549373 λ=-0.1 + 0.0im
iter 2 err:0.8397508140476416 λ=-0.6418389474323298 + 0.0im
iter 3 err:0.17336372619725743 λ=-0.08765753239354723 + 0.0im
iter 4 err:0.0005771170619943501 λ=-0.1029135620110966 + 0.0im
iter 5 err:4.762006833879597e-7 λ=-0.10285411985934721 + 0.0im
iter 6 err:4.074039107701665e-7 λ=-0.10285421074175707 + 0.0im
iter 7 err:2.6448037288912206e-8 λ=-0.10285417155884034 + 0.0im
iter 8 err:1.3926542408883378e-9 λ=-0.10285416898178967 + 0.0im
iter 9 err:6.324560618281378e-11 λ=-0.10285416884505445 + 0.0im
```

### User defined error 2: A user defined type

Due to the multiple dispatch and handling of types in Julia, code
may run faster if one uses types instead of function handles. It is
possible to do the same simulation as above with a user defined
type.

You first need to define a new type
```julia-repl
julia> struct RefErrmeasure <: Errmeasure; end
```
The error measure should then provided in the function
`estimate_error` which we now define as the relative eigenvalue error:
```julia-repl
julia> v0=ones(size(nep,1));
julia> (λref,v)=resinv(nep,v=v0,λ=230^2+1im,logger=0);
julia> function NonlinearEigenproblems.estimate_error(e::RefErrmeasure,λ,v)
         return abs(λ-λref)/abs(λ);
       end
julia> (λ,v)=resinv(nep,v=v0,λ=250^2+1im,logger=1,tol=1e-10,errmeasure=RefErrmeasure());
iter 1 err:0.12740916184575013 λ=62500.0 + 1.0im
iter 2 err:0.9535794095609479 λ=1.175146205389422e6 + 6663.738915456258im
iter 3 err:2.4041334742321228 λ=-38815.06145769358 + 1301.0064120392753im
....
iter 50 err:7.608430406801785e-11 λ=54550.13915567685 + 459.51715820395947im
```
The printouts of the last call correspond to the relative eigenvalue error.


## As a NEP-solver developer

NEP-solvers should use the `Errmeasure` as follows. The NEP-solver should take
as input an object of the type `Errmeasure`  or function. The fact
that it can be different types, is transparent and a NEP-solver
developer does not have to do anything to take care of that
if the following procedure is followed.

Suppose your solver is defined in a function with this
signature:
```julia
function mysolver(nep::NEP;errmeasure::ErrmeasureType=DefaultErrmeasure(nep))
```

In the main for loop you want to call the `estimate_error` function:

```julia
for k=1:maxit
    err=estimate_error(errmeasure,λ,v)
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
DefaultErrmeasure
```
```@docs
StandardSPMFErrmeasure
```
```@docs
ResidualErrmeasure
```
```@docs
EigvalReferenceErrmeasure
```
```@docs
estimate_error
```
```@docs
ErrmeasureType
```
