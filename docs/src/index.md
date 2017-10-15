
# NEPPACK 

Welcome to the NEPPACK documentation pages.

# Getting started


Let's solve a polynomial eigenvalue problem. First we need
to load the appropriate packages. 

```julia-repl
julia> using NEPCore, NEPSolver
```
The first time you run this, you will normally get an error message,
```julia-repl

```
NEPPACK builds on the community contributed julia
packages. These need to be installed the first time you use them,
by running for the corresponding packages:
```julia-repl
julia> Pkg.add("IterativeSolvers")
```

We are now ready to create and solve a nonlinear eigenvalue problem:
```julia-repl
julia> A0=[1 3; 5 6]; A1=[3 4; 6 6]
julia> nep=PEP([A0,A1,eye(2)])
NEPTypes.PEP(2, Array{Float64,2}[[1.0 3.0; 5.0 6.0], [3.0 4.0; 6.0 6.0], [1.0 0.0; 0.0 1.0]])
julia> λ,v=polyeig(nep)
(Complex{Float64}[1.36267+0.0im, -0.824084+0.280682im, -0.824084-0.280682im, -8.7145+0.0im], Complex{Float64}[-1.0+0.0im 0.739183-0.196401im 0.739183+0.196401im 0.627138+0.0im; 0.821812+0.0im -0.501408-0.375337im -0.501408+0.375337im 1.0+0.0im])
```
You have now solved your first nonlinear eigenvalue
problem with NEPPACK.

If we have a solution, then M(λ) should be singular,
with a singular vector v such that M(λ)v=0:
```julia-repl
julia> λ1=λ[1]; v1=v[:,1];
julia> norm((A0+A1*λ1+eye(2)*λ1^2)*v1)/norm(v1)
1.1502634749464687e-14
```




Reproduce an example in DDE-BIFTOOL


# NEPCore

The basic class is the abstract class `NEP` which represents
a NEP. 

```@docs
NEP
```

## Accessing the NEP

The nonlinear eigenvalue problem is defined by the data
stored in the NEP-class, and the NEP-solvers access
the data mainly through three main functions, `compute_Mder`
`compute_Mlincomb` and `compute_MM`.


```@docs
compute_Mder(nep::NEPCore.NEP, λ::Number, i::Integer) 
compute_Mlincomb(nep::NEPCore.NEP, λ::Number, V; a)
compute_MM
```



# NEPTypes

In order to use the methods,
the user has the possibility to implement their own
problem-specific functions above, or use one of the predefined
types. The most common one is the `SPMF_NEP`.


```@docs
SPMF_NEP
```

In order to construct an `SPMF_NEP`, we need to provide
the matrices and the functions.

```@docs
SPMF_NEP(AA::Array, fii::Array{Function,1}) 
```

```@docs
compute_Mder
```


```@docs
compute_Mlincomb
```

Let's try an equation $x=f(x)$. \(a \ne 0\)

```math
x=x_1+1
```







# Gallery


A large number of examples are provided in the `nep_gallery`.

```julia-repl
julia> using Gallery
julia> nep=nep_gallery("dep0")
julia> λ,v=newton(nep)
(-0.3587189459686265 + 0.0im, Complex{Float64}[0.284742+0.0im, -0.143316+0.0im, 0.278378+0.0im, -0.5009+0.0im, -0.613634+0.0im])
julia> norm(compute_Mlincomb(nep,λ,v))
4.718447854656915e-16
```

If MATLAB and the Berlin-Manchester collection areinstalled,
we can access them with the GalleryNLEVP
(which does MATLAB-access through julia's MATLAB-package).

```julia-repl
julia> using GalleryNLEVP
julia> nep=nep_gallery(NLEVP_NEP,"hadeler")
julia> λ,v=quasinewton(nep,λ=0.2,displaylevel=1,maxit=20,tol=1e-10);
julia> norm(compute_Mlincomb(nep,λ,v))/norm(v)
9.698206079849311e-11
```

Problems loaded from the Berlin-Manchester collection are NEP-objects
where every call to access a function generates a call to an
underlying MATLAB-session. Some problems in the Berlin-Manchester collection
have native support in NEPPACK, i.e., avoiding a MATLAB-access in every call.
The native equivalent object is generated with  `nlevp_make_native`:
```juli-repl
julia> using GalleryNLEVP
julia> nep1=nep_gallery(NLEVP_NEP,"gun")
julia> nep2=nlevp_make_native(nep1);
julia> norm(compute_Mder(nep1,0)-compute_Mder(nep2,0),1)
0.0
```

Stand-alone implementation can be accessed in a similar way, e.g.,
a native implementation of the waveguide eigenvalue problem:
```julia-repl
julia> using GalleryWaveguide
julia> nep=nep_gallery(WEP,benchmark_problem="TAUSCH");
```





# Compiling the documentation

Compile this documentation page by running:
```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ julia --color=yes make.jl &&  mkdocs build --clean
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ firefox site/index.html
```
If you want this to appear on our documentation page
[https://gitr.sys.kth.se/pages/nep-pack/nep-pack-alpha/](https://gitr.sys.kth.se/pages/nep-pack/nep-pack-alpha/)
you need to push it to the `gh-branch`, e.g.,  by running
```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ export DOCSDIR=`pwd`
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ cd /tmp
jarl@bjork:/tmp$ git clone -b "gh-pages" git@gitr.sys.kth.se:nep-pack/nep-pack-alpha.git
jarl@bjork:/tmp$ cd nep-pack-alpha
jarl@bjork:/tmp/nep-pack-alpha$ cp -r $DOCSDIR/site/* .
jarl@bjork:/tmp/nep-pack-alpha$ git add *
jarl@bjork:/tmp/nep-pack-alpha$ git commit . -m "refresh docs"
jarl@bjork:/tmp/nep-pack-alpha$ git push
```


More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)




