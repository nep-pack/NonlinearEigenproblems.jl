
<a id='NEPPACK-1'></a>

# NEPPACK


NEPPACK is a package with implementations of methods to solve nonlinear eigenvalue problems of the type: Find $(λ,v)\in\mathbb{C}\times\mathbb{C}^n$ such that


$$
M(λ)v=0
$$


and $v\neq 0$. 


<a id='Getting-started-1'></a>

# Getting started


Download NEP-PACK and get into the correct directory and start julia


```
jarl@bjork:~/src/nep-pack-alpha$
jarl@bjork:~/src/nep-pack-alpha$ julia
   _       _ _(_)_     |  A fresh approach to technical computing
  (_)     | (_) (_)    |  Documentation: https://docs.julialang.org
   _ _   _| |_  __ _   |  Type "?help" for help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 0.6.2 (2017-12-13 18:08 UTC)
 _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release
|__/                   |  x86_64-pc-linux-gnu

julia> 
```


First we need to make julia find the NEP-PACK-files in its path


```julia-repl
julia> push!(LOAD_PATH, string(ENV["HOME"],"/src/nep-pack-alpha/src/"))
```


and then we can start to load the appropriate packages. 


```julia-repl
julia> using NEPSolver, NEPTypes
```


The first time you run this, you will normally get an error message,


```julia
ERROR: LoadError: LoadError: ArgumentError: Module IterativeSolvers not found in current path.
Run `Pkg.add("IterativeSolvers")` to install the IterativeSolvers package.
```


NEPPACK builds on contributed julia packages. These need to be installed the first time you use NEPPACK, by running for the corresponding packages:


```julia-repl
julia> Pkg.add("IterativeSolvers")
```


We are now ready to create and solve a nonlinear eigenvalue problem, in this illustrative example we set


$$
M(λ)=\begin{bmatrix}1&3\newline5&6\end{bmatrix}+
λ\begin{bmatrix}3&4\newline6&6\end{bmatrix}+
λ^2\begin{bmatrix}1&0\newline0&1\end{bmatrix}
$$


The following code creates this NEP (which is a so called polynomial eigenvalue problem) and solves it using the NEP solution method implemented in `polyeig()`:


```julia-repl
julia> A0=[1 3; 5 6]; A1=[3 4; 6 6]
julia> nep=PEP([A0,A1,eye(2)])
NEPTypes.PEP(2, Array{Float64,2}[[1.0 3.0; 5.0 6.0], [3.0 4.0; 6.0 6.0], [1.0 0.0; 0.0 1.0]])
julia> λ,v=polyeig(nep)
(Complex{Float64}[1.36267+0.0im, -0.824084+0.280682im, -0.824084-0.280682im, -8.7145+0.0im], Complex{Float64}[-1.0+0.0im 0.739183-0.196401im 0.739183+0.196401im 0.627138+0.0im; 0.821812+0.0im -0.501408-0.375337im -0.501408+0.375337im 1.0+0.0im])
```


You have now solved your first nonlinear eigenvalue problem with NEPPACK. 


In order to verify that we have a solution, we can check that  $M(λ)$ is singular, with a singular vector $v$ such that $M(λ)v=0$:


```julia-repl
julia> λ1=λ[1]; v1=v[:,1];
julia> norm((A0+A1*λ1+eye(2)*λ1^2)*v1)/norm(v1)
1.1502634749464687e-14
```


<a id='Accessing-more-complicated-applications-1'></a>

# Accessing more complicated applications


We have made benchmark examples available in the module `Gallery`. Use it by loading the module and calling the function `nep_gallery`:


```julia-repl
julia> using Gallery
julia> nep=nep_gallery("dep0",100);
julia> size(nep)
(100, 100)
julia> λ,v=mslp(nep,tol=1e-10);
julia> λ
0.23169217667341738 - 2.1866254654451488e-16im
julia> size(v)
(100,)
julia> resnorm=norm(compute_Mlincomb(nep,λ,v))
3.124042808475689e-14
```


Information about the gallery can be found by typing `?nep_gallery`. The second arument in `nep_gallery` is a problem parameter, in this case specifying that the  size of the problem should be `100`. The example solves the problem with the method MSLP. The parameter `tol` specifies the tolerance for iteration termination. Type `?mslp` for more information about this method.


Under construction: Reproduce an example in DDE-BIFTOOL. Here is a Benchmark Example. 


<a id='Compiling-the-documentation-1'></a>

# Compiling the documentation


Compile this documentation page by running:


```
jarl@bjork:~/src/nep-pack-alpha/docs$ julia --color=yes make.jl &&  mkdocs build --clean
jarl@bjork:~/src/nep-pack-alpha/docs$ firefox site/index.html
```


If you want this to appear on our documentation page [https://nep-pack.github.io/nep-pack-alpha/](https://nep-pack.github.io/nep-pack-alpha/) you need to push it to the `gh-branch`, e.g.,  by running


```
jarl@bjork:~/src/nep-pack-alpha/docs$ export DOCSDIR=`pwd`
jarl@bjork:~/src/nep-pack-alpha/docs$ cd /tmp
jarl@bjork:/tmp$ git clone -b "gh-pages" git@github.com:nep-pack/nep-pack-alpha.git
jarl@bjork:/tmp$ cd nep-pack-alpha
jarl@bjork:/tmp/nep-pack-alpha$ cp -r $DOCSDIR/site/* .
jarl@bjork:/tmp/nep-pack-alpha$ git add *;  git commit . -m "refresh docs"; git push
```


More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)

