# Logger


## Basic usage
NEP-PACK provides considerable functionality to control the printouts and information of the NEP-solvers.
All NEP-solvers take the keyword argument `logger` which specifies if things should be stored in
a logger and/or printed. The main loggers are the [`PrintLogger`](@ref) which
only provides printouts, and [`ErrorLogger`](@ref) which stores error information.

We illustrate with a combination with the [error measure](errmeasure.md). This example
shows how to plot the eigenvalue error of a NEP-solver by using a reference solution.

First let us only user `logger=1` in combination with a [`EigvalReferenceErrmeasure`](@ref).
```julia-repl
julia> A0=[3.0 4 ; 5 6]; A1=[-1.0 0 ; 3.0 1.0];
julia> nep=DEP([A0,A1]); # Delay eigenvalue problem
julia> (λref,_)=resinv(nep,v=[1;1],λ=8,tol=1e-16); # Compute a reference solution
julia> resinv(nep,v=[1;1],λ=8,logger=1,errmeasure=EigvalReferenceErrmeasure(nep,λref));
Precomputing linsolver
iter 1 err:1.2171484853011378 λ=8.0 + 0.0im
iter 2 err:0.21696340485295096 λ=9.000185080448187 + 0.0im
iter 3 err:0.032989925133875886 λ=9.250138410435014 + 0.0im
iter 4 err:0.004864643426348181 λ=9.21228384187479 + 0.0im
iter 5 err:0.0007206309370317854 λ=9.21786911623817 + 0.0im
iter 6 err:0.00010667933045382938 λ=9.217041805970684 + 0.0im
iter 7 err:1.579396864670457e-5 λ=9.217164279269785 + 0.0im
iter 8 err:2.3382761789036977e-6 λ=9.217146147024959 + 0.0im
iter 9 err:3.461794584325162e-7 λ=9.217148831480596 + 0.0im
iter 10 err:5.1251506150151727e-8 λ=9.217148434049632 + 0.0im
iter 11 err:7.587733108493921e-9 λ=9.217148492888871 + 0.0im
iter 12 err:1.1233556307388426e-9 λ=9.217148484177782 + 0.0im
iter 13 err:1.6631140908884845e-10 λ=9.21714848546745 + 0.0im
iter 14 err:2.4622082150926872e-11 λ=9.217148485276516 + 0.0im
iter 15 err:3.645084234449314e-12 λ=9.217148485304783 + 0.0im
iter 16 err:5.400124791776761e-13 λ=9.217148485300598 + 0.0im
iter 17 err:7.993605777301127e-14 λ=9.217148485301218 + 0.0im
iter 18 err:1.0658141036401503e-14 λ=9.217148485301127 + 0.0im
```
The displayed err are eigenvalue errors and we now wish to plot them:
```jula-repl
julia> logger=ErrorLogger();
julia> resinv(nep,v=[1;1],λ=8,
    errmeasure=EigvalReferenceErrmeasure(nep,λref),logger=logger);
julia> errvec=logger.errs[1:17,1]; # This contains the iteration error
```
We use `Plots` for plotting:
```julia-repl
julia> using Plots;
julia> plot(errvec,yaxis=:log,marker=:star,xlabel="iteration",ylabel="eigval error")
```
The theory predicts linear convergence, which we also observe.
```@raw html
<br>
<img src="https://nep-pack.github.io/NonlinearEigenproblems.jl/logger_resinv_conv.png" height=300>
```

## Logger types
```@docs
Logger
PrintLogger
ErrorLogger
```

## Advanced usage

The logging functionality can be extended in case you want to collect
(or throw away) some of the information.
You need to create a new type which implements the following methods.

```@docs
push_info!
push_iteration_info!
```
