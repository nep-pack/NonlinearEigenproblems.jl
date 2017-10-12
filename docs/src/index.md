
# NEPPACK documentation

This is the documentation of the NEPPACK package.



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



# NEP methods

## Newton type methods
```@docs
NEPSolver.newton
NEPSolver.augnewton
NEPSolver.resinv
NEPSolver.quasinewton
NEPSolver.mslp
NEPSolver.rfi
NEPSolver.newtonqr
NEPSolver.implicitdet
```
## Projection methods
```@docs
NEPSolver.nlar
```
## Arnoldi type methods
```@docs
NEPSolver.iar
NEPSolver.tiar
NEPSolver.infbilanczos
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



