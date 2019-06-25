
# The Gallery function
A large number of examples are provided in the `nep_gallery`.

```julia-repl
julia> using Gallery
julia> nep=nep_gallery("dep0")
julia> λ,v=newton(nep)
(-0.3587189459686265 + 0.0im, Complex{Float64}[0.284742+0.0im, -0.143316+0.0im, 0.278378+0.0im, -0.5009+0.0im, -0.613634+0.0im])
julia> norm(compute_Mlincomb(nep,λ,v))
4.718447854656915e-16
```

```@docs
nep_gallery
```

# Berlin-Manchester collection
If MATLAB and the [Berlin-Manchester collection](http://www.maths.manchester.ac.uk/our-research/research-groups/numerical-analysis-and-scientific-computing/numerical-analysis/software/nlevp/) are installed,
we can access them with the GalleryNLEVP
(which does MATLAB-access through Julia's MATLAB-package).

```julia-repl
julia> using GalleryNLEVP
julia> nep=nep_gallery(NLEVP_NEP,"hadeler")
julia> λ,v=quasinewton(nep,λ=0.2,logger=1,maxit=20,tol=1e-10);
julia> norm(compute_Mlincomb(nep,λ,v))/norm(v)
9.698206079849311e-11
```

Problems loaded from the Berlin-Manchester collection are NEP-objects
where every call to access a function generates a call to an
underlying MATLAB-session. Some problems in the Berlin-Manchester collection
have native support in NEP-PACK, i.e., avoiding a MATLAB-access in every call;
see `nep_gallery` above.

# Other gallery examples
Stand-alone implementation can be accessed in a similar way, e.g.,
a native implementation of the Waveguide Eigenvalue Problem can be accessed as
```julia-repl
julia> using GalleryWaveguide
julia> nep=nep_gallery(WEP,benchmark_problem="TAUSCH");
```
