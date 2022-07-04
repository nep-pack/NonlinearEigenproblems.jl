# Tutorial: Implementing your own method

Although we try to provide state-of-the-art algorithms
in NEP-PACK, you may want to implement a solver
which is not available in NEP-PACK.
By using the NEP-PACK data types and structures when you implement your solver,
you can make your life easier in several ways.
You do not need to know the internals of NEP-PACK.
Correct usage will give you access to many applications,
helper functionality to combine with,
and you will have to possibility to compare your method
with other solvers.
We now illustrate how to implement your
own NEP-solver.

## Halley's method

[Halley's method for root-finding of nonlinear scalar equations](https://en.wikipedia.org/wiki/Halley%27s_method)
has fast local convergence - even faster
than Newton's method in terms convergence order, and
often faster in terms of number of iterations.
A NEP can be formulated as a
root-finding problem since a solution will always
satisfy
```math
f(λ)=\det(M(λ))=0
```
The application of Halley's method to this nonlinear scalar equation
will serve as an example solver, although it does, to our
knowledge, not lead to a competitive algorithm.
Halley's method for the root-finding problem is
defined by the iteration
```math
λ_{k+1}=λ_k-\frac{2f(λ_k)f'(λ_k)}{2(f'(λ_k))^2-f(λ_k)f''(λ_k)}
```
There are formulas for the
derivatives of the determinant, we will here for
simplicity just use finite difference approximation to
estimate the derivatives, i.e.,
```math
 f'(λ)\approx \frac{f(λ+δ)-f(λ-δ)}{2δ}
```
```math
 f''(λ)\approx \frac{f(λ+δ)-2f(λ)+f(λ-δ)}{δ^2}
```
## Implementation in NEP-PACK (preliminary version)

Let us first define our solver function
and introduce the function whose roots we wish to find.
The matrix ``M(λ)`` is obtained by a call to the
[`compute_Mder`](@ref)-function.
```julia
using NonlinearEigenproblems
function halley(nep::NEP;λ=0.0,δ=sqrt(eps()),maxit=100,tol=eps()*100)
   f=s-> det(compute_Mder(nep,s)); # The objective function
   # More code here
end
```
The main loop (which should go in `# More code here`) can be implemented,
in a way that does not involve many function
evaluations, as follows:
```julia
   for i=1:maxit
       fλ=f(λ)
       fλp=f(λ+δ)
       fλm=f(λ-δ)
       fp=(fλp-fλm)/(2δ)
       fpp=(fλp-2*fλ+fλm)/(δ^2)
       Δλ=2*fλ*fp/(2*fp^2-fλ*fpp);
       λ=λ-Δλ;
       @show (i,λ)
       if (abs(Δλ)<tol)
          return λ
       end
   end
```
Let us now test the code on a benchmark problem:
```julia
julia> nep=nep_gallery("dep0");
julia> λ=halley(nep)
(i, λ) = (1, -0.13876571372157542)
(i, λ) = (2,  0.15938372556136426)
(i, λ) = (3, -0.15955391446207692)
(i, λ) = (4, -0.15955391823299248)
```
Clearly, the algorithm terminates after 4 iterations.
We can verify that this is actually
a solution easily if we also
have an approximate eigenvector. An eigenvector
can be computed/estimated by essentially one step of inverse iteration,
on the matrix ``M(λ)``:
```julia
julia> x=normalize(compute_Mder(nep,λ)\ones(size(nep,1)))
5-element Array{Float64,1}:
  0.14358324743994907
  0.9731847884093298
 -0.12527093093249475
  0.031821422867456914
  0.12485915894832478
```
The residual norm  ``||M(λ)x||`` does indeed become almost zero
so it seems we have a solution:
```julia
julia> norm(compute_Mlincomb(nep,λ,x))
7.093661646042283e-16
```

## Implementation in NEP-PACK (full version)

In the following we illustrate a more advanced
usage of the NEP-PACK method development:
NEP-PACKs logging facility  and error estimation.
See [`Logger`](logger.md) and [`Errmeasure`](errmeasure.md). This
gives access
to other ways to measure error as well as a logging and
inspection of error history in a way that is
the same for all solvers and simplifies
comparisons.

```julia
using NonlinearEigenproblems, LinearAlgebra, Plots
function halley(nep::NEP;λ=0.0,δ=sqrt(eps()),maxit=100,
                tol=eps()*100,logger=0,
                errmeasure = DefaultErrmeasure(nep))
    # Setup the logger.
    @parse_logger_param!(logger);

    n=size(nep,1);
    f=s-> det(compute_Mder(nep,s)); # The objective function


    for i=1:maxit
        fλ=f(λ)
        fλp=f(λ+δ)
        fλm=f(λ-δ)
        fp=(fλp-fλm)/(2δ)
        fpp=(fλp-2*fλ+fλm)/(δ^2)
        Δλ=2*fλ*fp/(2*fp^2-fλ*fpp);
        λ=λ-Δλ;
        # Compute an eigenvector. This will not work if the
        # eigenvector is orthogonal to ones(n)
        x=normalize(compute_Mder(nep,λ)\ones(n));
        err=estimate_error(errmeasure,λ,x)  # Estimate the error
        push_iteration_info!(logger,i; λ=λ,err=err) # Put it into the log
        if (err<tol)
            return (λ,x)
        end
    end
end
```

We can now run our new method using
with a `logger=1` keyword argument
so we get the standardized output of iteration info:
```julia-repl
julia> (λ,x)=halley(nep,logger=1);
iter 1 err:0.010384216303530201 λ=-0.13876571372157542
iter 2 err:8.082978338039669e-5 λ=-0.15938372556136426
iter 3 err:1.7901681647471861e-9 λ=-0.15955391446207692
iter 4 err:1.0389976569127096e-16 λ=-0.15955391823299248
julia> norm(compute_Mlincomb(nep,λ,x))
7.093661646042283e-16
```
The use of the NEP-PACK logging functionality makes it
very easy to visualize the error. If you now want to plot the error history,
you can use the [`ErrorLogger`](@ref):
```julia-repl
julia> mylogger=ErrorLogger()
julia> (λ,x)=halley(nep,logger=mylogger);
julia> plot(mylogger.errs[1:10,1],yaxis=:log)
```
We clearly observe the superlinear convergence:
```@example
using Gadfly # hide
z=[ 0.08492602120772309   # hide
        0.07450867012944977 # hide
        0.032639292900081246 # hide
        0.00281602165251169 # hide
        1.1025990567599428e-5 # hide
        1.0638098128402615e-10 # hide
        4.942402279980973e-17 # hide
       ]; # hide
plot(y=z, Scale.y_log10()) # hide
```

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_NEWMETHOD)
