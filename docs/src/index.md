
# NEPPACK 

NEPPACK is a package with implementations of methods to solve nonlinear eigenvalue problems of
the type: Find ``(λ,v)\in\mathbb{C}\times\mathbb{C}^n`` such that
```math
M(λ)v=0
```
and ``v\neq 0``. 


# Getting started


First we need to load the appropriate packages. 
```julia-repl
julia> using NEPCore, NEPSolver
```
The first time you run this, you will normally get an error message,
```julia
ERROR: LoadError: LoadError: ArgumentError: Module IterativeSolvers not found in current path.
Run `Pkg.add("IterativeSolvers")` to install the IterativeSolvers package.
```
NEPPACK builds on contributed julia packages.
These need to be installed the first time you use them,
by running for the corresponding packages:
```julia-repl
julia> Pkg.add("IterativeSolvers")
```

We are now ready to create and solve a nonlinear eigenvalue problem, in this example

```math
M(λ)=\begin{block}1&3\\5&6\end{block}+
λ\begin{block}3&4\\6&6\end{block}+
λ^2\begin{block}1&0\\0&1\end{block}
```
The following code creates this NEP (which is a so called polynomial eigenvalue problem)
and solves it using the NEP solution method implemented in `polyeig()`:
```julia-repl
julia> A0=[1 3; 5 6]; A1=[3 4; 6 6]
julia> nep=PEP([A0,A1,eye(2)])
NEPTypes.PEP(2, Array{Float64,2}[[1.0 3.0; 5.0 6.0], [3.0 4.0; 6.0 6.0], [1.0 0.0; 0.0 1.0]])
julia> λ,v=polyeig(nep)
(Complex{Float64}[1.36267+0.0im, -0.824084+0.280682im, -0.824084-0.280682im, -8.7145+0.0im], Complex{Float64}[-1.0+0.0im 0.739183-0.196401im 0.739183+0.196401im 0.627138+0.0im; 0.821812+0.0im -0.501408-0.375337im -0.501408+0.375337im 1.0+0.0im])
```
You have now solved your first nonlinear eigenvalue
problem with NEPPACK. 

In order to verify that we have a solution, we can check that  ``M(λ)`` is singular,
with a singular vector ``v`` such that ``M(λ)v=0``:
```julia-repl
julia> λ1=λ[1]; v1=v[:,1];
julia> norm((A0+A1*λ1+eye(2)*λ1^2)*v1)/norm(v1)
1.1502634749464687e-14
```




Reproduce an example in DDE-BIFTOOL

Here is a Benchmark Example..





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




