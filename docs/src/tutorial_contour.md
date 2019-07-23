# Contour integral tutorial


NEP-PACK contains several implementations of different
contour integral methods.
Although they have been worked out and
presented independently (in different articles),
we have implemented them in a way where the
common componentsare the same.

Contour integral methods are have one property
which makes them attractive from the perspective of
parallelization, which we will illustrate in
the final example below.


## Basic usage

We illustrate both how to make
call to contour Beyn's contour
integral method and the block SS method
of Asakura and Sakurai.

First set up a large and sparse problem:
```julia
julia> using SparseArrays;
julia> n=100;
julia> A0=spdiagm(0 => ones(n))
julia> A1=spdiagm(-2 => ones(n-2), 1 => (1:n-1)/n)
julia> A2=spdiagm(-1 => ones(n-1), 0 => (1:n)/n, 1 => sin.(range(0,10,length=n-1)))
julia> nep=SPMF_NEP([A0,A1,A2],[s->one(s), s->s^2, s->sin(s)])
```
and call the two integral solution methods:
```julia
julia> (λ,v)= contour_beyn(nep,radius=0.5,k=10);
```
and we can verify that we found some good solutions
```julia
julia> λ
2-element Array{Complex{Float64},1}:
  -0.2584032183498001 + 0.41445913392459055im
 -0.25840321834980007 - 0.4144591339245906im
julia> norm(compute_Mlincomb(nep,λ[1],normalize(v[:,1])))
2.207347574532086e-15
julia> norm(compute_Mlincomb(nep,λ[2],normalize(v[:,2])))
2.052254464542942e-15
```
For later comparison we also use `contour_block_SS`
```julia
julia> (λ,v)= contour_block_SS(nep,radius=0.5,k=10);
julia> λ
4-element Array{Complex{Float64},1}:
  -0.2584032183498005 + 0.4144591339245909im
  -0.3635139038679691 + 0.33597403312059815im
 -0.36351390386796917 - 0.335974033120598im
 -0.25840321834980057 - 0.4144591339245917im
```
and the corresponding residual norms
```julia
julia> norm(compute_Mlincomb(nep,λ[1],normalize(v[:,1])))
2.9727040010325543e-15
julia> norm(compute_Mlincomb(nep,λ[2],normalize(v[:,2])))
2.5579929787574954e-15
julia> norm(compute_Mlincomb(nep,λ[3],normalize(v[:,3])))
2.360173099697166e-15
julia> norm(compute_Mlincomb(nep,λ[4],normalize(v[:,4])))
3.920633762821404e-15
```
The functions `contour_beyn` and `contour_block_SS`
have compatible keyword argumengs. The value `radius=0.5`,
means that we numerically integrate  a circle of radius `0.5`,
centered at `σ=0`. We should expect the method to find
eigenvalues (hopefully all eigenvalues) within that disk.
(Our implementation also supports ellipses, by specifying
`radius` as a length two vector with the two radii of the ellipse.)
The value `k=5` specifies how many columns the search subspace has.
In general, we do not obtain more `k` eigenvalues.

It seems that in this case `contour_block_SS` is better
since it finds two eigenvalues inside the disk which
`contour_beyn` misses. However, a closer look reveals
that the call to  `contour_block_SS` also requires
more computation time.

## Your own integration method

The contour integral methods are based on numerical quadrature.
There are many different ways to carry out quadrature,
and NEP-PACK provides a way to use user-defined
quadrature methods.
The default is to use the trapezoidal rule. When we parameterize
a circle (or ellipse) with a phase, the integrand is periodic
and the trapezoidal rule works particularly well.
It is however not the only option for quadrature and
we can for instance implement a gauss quadrature,
in this case by using the functionality in the package `QuadGK`:
```julia
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

```julia
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
as a parameter in the call after the `nep`:
```julia
julia> (λ,v)= contour_block_SS(nep,GaussIntegrator,radius=0.5, k=10);
julia> λ
4-element Array{Complex{Float64},1}:
 -0.25840321834980035 + 0.414459133924587im
    -0.36351390386797 + 0.3359740331206035im
  -0.3635139038679633 - 0.33597403312059543im
  -0.2584032183498066 - 0.4144591339245837im
```
Let's make it print some pretty decoration
during the progress of the method.
In the code where it currently says
`## Extra code goes here` we will now insert
```julia
if (mod(i,round(N/50))==1)
   print(".")
end
``
and `println()` in the second code insertion.
This way, we will print a progress bar, which
prints in total (approximately) 50 dots.
You will see dots gradually appearing:
```julia
julia> (λ,v)= contour_beyn(nep,GaussIntegrator,radius=0.5,k=10);
..................................................
```


## Parallellized integration method

The main computational effort of the contour
integral methods lies in solving many linear systems.
This is done in the call to `f` in
the `integrate_interval`-function. Since they are completely
independent operation in the for-loop, they can
be easily parallelized.