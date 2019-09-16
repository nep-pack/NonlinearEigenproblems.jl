# Projection

Many NEP-solvers are based on a computation of a solution
to a projected problem, i.e., if ``V,W\in\mathbb{R}^{n\times p}``
we need to solve the (smaller) NEP
```math
W^HM(λ)Vz=0
```
This is sometimes called a nonlinear Rayleigh-Ritz procedure,
or a direct projection. These are *inner solvers* for many NEP-solvers.


NEP-PACK provides a framework to handle projected problems
and inner solves. This is implemented
into two separate components:

* [Projection](@ref): As a user (or NEP-solver developer) you can create a new object corresponding to the projection. In NEP-PACK, the projection is again an object of with type inheriting from [`NEP`](@ref). More precisely, it is a [`Proj_NEP`](@ref) which you normally create with the function [`create_proj_NEP`](@ref).
* [Inner solvers](@ref): Since the projected problem is again a `NEP`, in principle any of the NEP-solvers of this package can be used. This is handled by the `InnerSolver` objects which are wrappers for corresponding NEP-solvers such that we can pass appropriate parameters to the inner soler. The inner solver is controlled by the `inner_solver_method` keyword in many NEP-solvers. By default [`DefaultInnerSolver`](@ref) is used.



As a NEP-user, you often do not need to care about how the
projection is handled, e.g., if you use the type [`SPMF_NEP`](@ref)
with only a few terms. For instance,
if you wish to use the infinite Arnoldi method ([`iar`](@ref))
to handle the project solves in the nonlinear
Arnoldi method ([`nlar`](@ref)), you can call `nlar` with the kwarg `inner_solver_method=`[`IARInnerSolver`](@ref)`()`:

```julia-repl
julia> nep=nep_gallery("dep0_tridiag");
julia> λ,v=nlar(nep,neigs=1,inner_solver_method=IARInnerSolver(),logger=1);
Using inner solver IARInnerSolver(1.0e-13, 80, :ones, false, NonlinearEigenproblems.NEPSolver.iar)
iter 1 err:0.07075594099201046 λ=-0.09271600160844096 + 0.03235352263155604im
iter 2 err:0.03648528716335285 λ=0.1565592513189133 + 0.013613366856692472im
iter 3 err:0.026208926915836202 λ=-0.03227579276149921 - 0.02208846862836847im
iter 4 err:0.007172980379743359 λ=-0.05583493095814305 - 0.003087442900200613im
iter 5 err:0.006924911617556792 λ=-0.05401508521132263 + 0.001517986336058301im
iter 6 err:0.002090878656094346 λ=-0.05536501763852141 + 0.00024284078743963842im
iter 7 err:0.0006103719429985026 λ=-0.05575237680802954 - 0.00013474280139834488im
iter 8 err:0.0002665363571995074 λ=-0.05558344084437247 + 6.558335069789856e-6im
iter 9 err:0.023093036823993347 λ=-0.02304829345220316 + 0.0005911146839749534im
iter 10 err:1.092985097378946e-5 λ=-0.055653357185449566 - 2.575336159079228e-6im
iter 11 err:5.316449394914625e-7 λ=-0.05565515933487699 + 2.446717728666673e-7im
iter 12 err:0.016299715925824968 λ=0.023211712543947764 - 0.033813269858071315im
iter 13 err:4.899222816433482e-8 λ=-0.05565520569884035 - 2.0175841868490324e-8im
iter 14 err:4.320444084682331e-9 λ=-0.055655212526421895 + 1.2690665519371834e-9im
iter 15 err:7.639277885601529e-10 λ=-0.05565521363070881 + 3.971813908747728e-11im
iter 16 err:4.941632484271334e-11 λ=-0.055655213920101095 - 7.493472629149413e-12im
iter 17 err:0.010663980713333146 λ=0.008967322670902417 + 0.016345149953039387im
iter 18 err:3.2317477531099844e-12 λ=-0.05565521393097803 - 6.246051566770699e-13im
iter 19 err:2.3964655361108506e-13 λ=-0.05565521393184677 + 6.221189257479347e-14im
iter 20 err:1.1735483724254833e-13 λ=-0.055655213931828935 + 2.058434802154811e-15im
iter 21 err:8.164760090914088e-15 λ=-0.055655213931847865 + 7.498992677576017e-16im
****** 1 converged to eigenvalue: -0.055655213931847865 + 7.498992677576017e-16im errmeasure:8.164760090914088e-15
```

The logging of the inner solver is controlled by the kwarg `inner_logger`,
which follows the same framework as the standard [NEP-PACK Logger](logger.md).
This produces very verbose output illustrating
also the convergence of the inner solve:
```julia-repl
julia> λ,v=nlar(nep,neigs=1,inner_solver_method=IARInnerSolver(),logger=1,inner_logger=1);

Using inner solver IARInnerSolver(1.0e-13, 80, :ones, false, NonlinearEigenproblems.NEPSolver.iar)
-
--
---
----
-----
------
=------
+-------
iter 1 err:0.06907648709827012 λ=-0.13302652304722704 + 0.0499583011875092im
-
--
---
----
-----
------
=------
+-------
iter 2 err:0.03702238592099922 λ=0.08696844790344242 + 0.010741310688729204im
-
--
---
----
-----
------
-------
=-------
+=-------
iter 3 err:0.029408773621139236 λ=0.0076466477038438325 - 0.07172981749577159im
-
--
---
----
-----
------
-------
--------
+--------
...
```

[Rayleigh functional computation](@ref), which corresponds to projection
with $p=1$, is also handled with this framework.

## Inner solvers

The inner solvers inherit from [`InnerSolver`](@ref).
The following inner solvers are available by default.


```@docs
NewtonInnerSolver
```

```@docs
IARInnerSolver
```

```@docs
IARChebInnerSolver
```

```@docs
ContourBeynInnerSolver
```

```@docs
PolyeigInnerSolver
```

```@docs
SGIterInnerSolver
```

```@docs
NleigsInnerSolver
```

```@docs
DefaultInnerSolver
```


## Inner solvers: Advanced usage

You can define your own inner solver by
inheriting from `InnerSolver` and implementing
the function `inner_solve`. Since the `inner_solve`
obtains information from the solver via
keyword arguments, you need to end your
method signature with `kwargs...)`.

```@docs
InnerSolver
```

```@docs
inner_solve
```


## Projection

The NEP-PACK functionality for projected problems
are represented by [projection types](@ref).
Normally, the projection is created by
[`create_proj_NEP`](@ref) from a standard NEP.
After creating a projected NEP, you can set
the projection subspace (represented by the
matrices `V` and `W`) using
[`set_projectmatrices!`](@ref) or
[`expand_projectmatrices!`](@ref).
```julia-repl
julia> A=[1 0 0; 0 1.0 0; 0 0 1]; B=[1 2 3; 3 3 3 ; 4 -1 -1.0];
julia> nep=SPMF_NEP([A, B], [s->s, s->s^5]);
julia> pnep=create_proj_NEP(nep);
julia> W=[4 1 ; 6 1  ; 6.0 2]; V=[3 3;3 4.0;4.0 -1];
julia> set_projectmatrices!(pnep,W,V); # modifies pnep
julia> λ=3.0+1im;
julia> W'*compute_Mder(nep,λ)*V
2×2 Array{Complex{Float64},2}:
 -3366.0+92958.0im  -2238.0+61334.0im
  -690.0+19290.0im   -513.0+13909.0im
julia> compute_Mder(pnep,λ)
2×2 Array{Complex{Float64},2}:
 -3366.0+92958.0im  -2238.0+61334.0im
  -690.0+19290.0im   -513.0+13909.0im
```
Effectively, the `Proj_NEP` creates [compute functions](compute_functions.md),
which are designed to be as efficient as possible.

### Projection functions

You can create a projected NEP with `create_proj_NEP`, and specify the
projection space with 
[`set_projectmatrices!`](@ref) and [`expand_projectmatrices!`](@ref).

```@docs
create_proj_NEP
```

```@docs
set_projectmatrices!
```

```@docs
expand_projectmatrices!
```

### Projection types
NEPs for which this projection can be computed
inherit from `ProjectableNEP`.

```@docs
ProjectableNEP
```

The result of the projection is represented in a `Proj_NEP`.

```@docs
Proj_NEP
```

One explicit instance is the `Proj_SPMF_NEP`.

```@docs
Proj_SPMF_NEP
```




## Rayleigh functional computation


```@docs
compute_rf
```

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_INNERSOLVE)
