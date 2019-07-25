# Contour integral tutorial


NEP-PACK contains several implementations of methods
in the family of approaches based on
contour integration.
Although they have been worked out and
presented independently
(in different research articles by different research groups),
we have implemented them in a unified
and extendible way.

Contour integral methods have one property
which makes them attractive from the perspective of
parallelization, which we will illustrate in
the final example below.


## Basic usage

The most popular methods contour
integral methods are Beyn's contour
integral method and the block SS method
of Asakura and Sakurai. We illustrate both of them.
First set up a large and sparse problem:
```julia-repl
julia> using SparseArrays;
julia> n=1000;
julia> A0=spdiagm(0 => ones(n))
julia> A1=spdiagm(-2 => ones(n-2), 0 => 30*(n:-1:1)/n,  1 => 3*ones(n-1))/3
julia> A2=spdiagm(-1 => ones(n-1), 0 => (1:n)/n, 1 => sin.(range(0,5,length=n-1)))/10
julia> nep=SPMF_NEP([A0,A1,A2],[s->one(s), s->s, s->exp(-s)])
```
and call the two integral solution methods:
```julia-repl
julia> (λ,v)= contour_beyn(nep,radius=0.5,k=10);
```
We can verify that we found some good solutions
```julia-repl
julia> λ
2-element Array{Complex{Float64},1}:
 -0.4938003805961036 + 0.03369433628038132im
 -0.4984653501095431 - 0.013414744968396205im
julia> norm(compute_Mlincomb(nep,λ[1],normalize(v[:,1])))
2.8693125572899838e-6
julia> norm(compute_Mlincomb(nep,λ[2],normalize(v[:,2])))
3.0028543096707394e-6
```
For comparison we also use `contour_block_SS`
```julia-repl
julia> (λ,v)= contour_block_SS(nep,radius=0.5,k=10);
julia> julia> λ
7-element Array{Complex{Float64},1}:
 -0.49789562317811836 + 0.029382625854591973im
  -0.5020899933123398 - 0.027308288264250524im
  -0.5006296180796399 + 0.011976675667372098im
  -0.5000287784310599 - 0.010301420892154335im
  -0.5044451294089868 - 0.0074606034247795975im
  -0.5001550771105308 - 0.00026147429323077303im
 -0.49957316937095864 + 0.003511006328045692im
```
and the corresponding residual norms
```julia-repl
julia> for j=1:7; @show norm(compute_Mlincomb(nep,λ[j],normalize(v[:,j]))); end
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 2.8693125572899838e-6
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 3.0028543096707394e-6
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 1.1514402700870265e-7
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 4.123810796391466e-8
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 6.261761794674978e-8
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 1.7269388863226059e-9
norm(compute_Mlincomb(nep, λ[j], normalize(v[:, j]))) = 2.994085882385125e-9
```
The functions `contour_beyn` and `contour_block_SS`
have compatible keyword argumengs. The kwarg `radius=0.5`,
means that we numerically integrate  a circle of radius `0.5`.
The center of the circle is given by the `σ`, argument and
by default `σ=0`. We should expect the method to find
eigenvalues (hopefully all eigenvalues) within that disk.
(Our implementation also supports ellipses, by specifying
`radius` as a length two vector with the two radii of the ellipse.)
The value `k=10` specifies how many columns the rectangular
probe matrix has.
In general, we do not obtain more `k` eigenvalues.

It seems that in this case `contour_block_SS` is better
since it finds eigenvalues  which
`contour_beyn` misses. However, a closer look reveals
that the additional eigenvalues
are outside the requested disc, and the
call to  `contour_block_SS` also requires
more computation time, making the comparison
unfair.

## Your own quadrature method

The contour integral methods are based on numerical quadrature.
There are many different ways to carry out quadrature,
and NEP-PACK provides a way to use user-defined
quadrature methods.
The default behaviour is to use the trapezoidal rule.
When we parameterize
a circle (or ellipse) with a phase, the integrand is periodic
and the trapezoidal rule works particularly well.
It is however not the only option for quadrature and
we can for instance implement a gauss quadrature,
in this case by using the functionality in the package `QuadGK`:
```julia-repl
julia> using Pkg
julia> Pkg.add("QuadGK");
julia> using QuadGK
```
The function `(x,w)=gauss(N)` provides weights and quadrature
points for a function to be integrated over the
interval `[-1,1]` with `N` quadrature points.

Before implementing the method, let us first have a look
at the documtation of `MatrixIntegrator`:

```@docs
MatrixIntegrator
```
Let us now combine the `gauss`-method in an implementation
of a numerical quadrature to be used in the quadrature
methods.

```julia-repl
julia> abstract type GaussIntegrator <: MatrixIntegrator; end
julia> import  NonlinearEigenproblems.NEPSolver.integrate_interval
julia> function integrate_interval(ST::Type{GaussIntegrator},::Type{T},f,gv,a,b,N,logger) where {T<:Number}
    x,w=gauss(N);        # Compute the Gauss weights
    w=w*(b-a)/2;         # Rescale w to interval [a,b]
    t=a .+ ((x .+ 1)/2)*(b-a); # Rescale t
    m=size(gv,1);
    # create the tensor and compute all quadratures
    S = zeros(T,size(f(t[1]))...,m)
    for i = 1:N
        ## Extra code goes here 1
        temp = f(t[i]) # Only computed once for all g-functions
        for j=1:m
            S[:,:,j] += temp*(gv[j](t[i])*w[i]);
        end
    end
    ## Extra code goes here 2
    return S
end
```
To specify this solver, you need to add the type you just created
as a parameter in the call. This is an argument (not a keyword argument)
after the argument `nep`:
```julia-repl
julia> (λ,v)= contour_block_SS(nep,GaussIntegrator,radius=0.5, k=10);
julia> λ
6-element Array{Complex{Float64},1}:
  -0.5030050924478993 + 0.025867789190345332im
  -0.4998917126923037 - 0.014647029189145597im
 -0.49991828738335686 - 0.007092586236661307im
  -0.5000067107140442 - 0.0026614262456865663im
 -0.49903549969757116 + 0.0075397370638041255im
   -0.501620024772268 + 0.00393810326235837im
```
Let's make it print some pretty decoration
during the progress of the method.
In the code where it currently says
`## Extra code goes here` we will now insert
```julia-repl
if (mod(i,round(N/50))==1)
   print(".")
end
```
and `println()` in the second code insertion.
This way, we will print a progress bar, which
prints in total (approximately) 50 dots.
You will see dots gradually appearing as a progress
bar:
```julia-repl
julia> (λ,v)= contour_beyn(nep,GaussIntegrator,radius=0.5,k=10);
..................................................
```


## Parallellized quadrature method

The main computational effort of the contour
integral methods lies in solving many linear systems.
This is done in the call to `f` in
the `integrate_interval`-function. Since they are completely
independent operations in the for-loop, they can
be easily parallelized.

Install the package `Distributed` and `BenchmarkTools` and include
with
```julia-repl
julia> using Distributed,BenchmarkTools
```
Similar to the previous example we make a new
type corresponding to our integrator and
explicitly import that
```julia-repl
julia> abstract type ParallelIntegrator <: MatrixIntegrator; end
julia> import  NonlinearEigenproblems.NEPSolver.integrate_interval
```
and define a function which computes the main for-loop in parallel using
the `@distributed` macro:
```julia-repl
julia> function integrate_interval(ST::Type{ParallelIntegrator},::Type{T},f,gv,a,b,N,logger) where {T<:Number}
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    m = size(gv,1);
    S = @distributed (+) for i = 1:N
        temp = f(t [i])
        Z=zeros(T,size(temp,1),size(temp,2),m);
        for j=1:m
            Z[:,:,j]=temp*gv[j](t [i]);
        end
        Z
    end
    return S
end
```
To use the parallelization you may need to start
julia with command-line arguments to specify the number
of parallel processes to be used, e.g., `-p 4`.
The `@btime` macro provides a way to measure how much faster
the parallel implementation is.
```julia-repl
julia> @btime (λ,v)= contour_block_SS(nep,ParallelIntegrator,radius=0.5, k=10);
  863.420 ms (1385 allocations: 10.46 MiB)
julia> @btime (λ,v)= contour_block_SS(nep,radius=0.5, k=10);
  2.990 s (321362 allocations: 5.84 GiB)
```
This is a speed up of 3.4, with `p=4` processes.