
<a id='NEPCore-Documentation-1'></a>

# NEPCore Documentation

<a id='NEPCore.NEP' href='#NEPCore.NEP'>#</a>
**`NEPCore.NEP`** &mdash; *Type*.



```
abstract NEP
```

A NEP object represents a nonlinear eigenvalue problem. All NEPs should implement

```julia-repl
size(nep::NEP,d)
```

and at least one of the following

```
M=compute_Mder(nep::NEP,λ::Number;der=0)
V=compute_Mlincomb(nep::NEP,λ::Number,V::Array{<:Number,2})
MM=compute_MM(nep::NEP,S::Array{<:Number,2},V::Array{<:Number,2})
```

<a id='NEPTypes.SPMF_NEP' href='#NEPTypes.SPMF_NEP'>#</a>
**`NEPTypes.SPMF_NEP`** &mdash; *Type*.



```
type SPMF_NEP <: AbstractSPMF
```

An SPMF_NEP is a NEP defined by a Sum of Products of Matrices and Functions, i.e.,

$$
M(λ)=∑_i A_i f_i(λ).
$$

All of the matrices $A_0,...$ are of size $n×n$ and $f_i$ are a functions. The  functions $f_i$ must be defined for matrices in the standard matrix function sense. 

<a id='NEPTypes.SPMF_NEP-Tuple{Array,Array{Function,1}}' href='#NEPTypes.SPMF_NEP-Tuple{Array,Array{Function,1}}'>#</a>
**`NEPTypes.SPMF_NEP`** &mdash; *Method*.



```
 SPMF_NEP(AA,fii,Schur_fact=false)
```

Creates a SPMF_NEP consisting of matrices `AA` and functions `fii`. `fii` must be an array of functions defined for matrices. `AA` is an array of matrices. `Schur_fact` specifies if the computation of `compute_MM` should be done by first pre-computing a Schur-factorization (which can be faster).

```julia-repl
julia> A0=[1 3; 4 5]; A1=[3 4; 5 6];
julia> id_op=S -> eye(S)
julia> exp_op=S -> expm(S)
julia> nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
julia> compute_Mder(nep,1)-(A0+A1*exp(1))
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0
```

<a id='NEPCore.compute_Mder' href='#NEPCore.compute_Mder'>#</a>
**`NEPCore.compute_Mder`** &mdash; *Function*.



```
compute_Mder(nep::NEP,λ::Number [,i::Integer=0])
```

Computes the ith derivative of `nep` evaluated in `λ`.

**Example**

This example shows that `compute_Mder(nep,λ,1)` gives the first derivative.

```julia-repl
julia> nep=nep_gallery("dep0");
julia> ϵ=1e-5;
julia> Aminus=compute_Mder(nep,λ-ϵ);
julia> Aminus=compute_Mder(nep,λ-ϵ);
julia> Aplus=compute_Mder(nep,λ+ϵ);
julia> norm((Aplus-Aminus)/(2ϵ)-compute_Mder(nep,λ,1))
1.990970375089371e-11
```

<a id='NEPCore.compute_Mlincomb' href='#NEPCore.compute_Mlincomb'>#</a>
**`NEPCore.compute_Mlincomb`** &mdash; *Function*.



```
compute_Mlincomb(nep::NEP,λ::Number,V;a=ones(size(V,2)))
```

Computes the linear combination of derivatives
$Σ_i a_i M^{(i)}(λ) v_i$

**Example**

This example shows that `compute_Mder` gives a result consistent with `compute_Mlincomb`. Note that `compute_Mlincomb` is in general faster since no matrix needs to be constructed.

```julia-repl
julia> nep=nep_gallery("dep0");
julia> v=ones(size(nep,1)); λ=-1+1im;
julia> norm(compute_Mder(nep,λ,1)*v-compute_Mlincomb(nep,λ,hcat(v,v),a=[0,1]))
1.0778315928076987e-15

```


```
compute_Mlincomb(nep::NEP,λ::Number,V,a::Array,startder::Integer)
```

Computes linear combination starting with derivative startder, i.e., $Σ_i a_i M^{(i+startder)}(λ) v_i$

The default implementation of this can be slow. Overload for specific NEP if you want efficiency (for aug_newton, IAR, ..).


Let's try an equation $x=f(x)$. (a \ne 0)


$$
x=x_1+1
$$


<a id='Compiling-the-documentation-1'></a>

# Compiling the documentation


Compile this documentation page by running:


```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ julia --color=yes make.jl
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ mkdocs build --clean
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ firefox site/index.html
```


If you want this to appear on our documentation page [https://gitr.sys.kth.se/pages/nep-pack/nep-pack-alpha/](https://gitr.sys.kth.se/pages/nep-pack/nep-pack-alpha/) you need to push it to the `gh-branch`, e.g.,  by running


```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ echo DOCSDIR=`pwd()`
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ cd /tmp
jarl@bjork:/tmp$ git clone -b "gh-pages" git@gitr.sys.kth.se:nep-pack/nep-pack-alpha.git
jarl@bjork:/tmp$ cd nep-pack-alpha
jarl@bjork:/tmp/nep-pack-alpha$ cp -r $DOCSDIR/* .
jarl@bjork:/tmp/nep-pack-alpha$ git add *
jarl@bjork:/tmp/nep-pack-alpha$ git commit .
jarl@bjork:/tmp/nep-pack-alpha$ git push
```


More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)


<a id='NEP-methods-1'></a>

# NEP methods


<a id='Newton-type-methods-1'></a>

## Newton type methods

<a id='NEPSolver.newton' href='#NEPSolver.newton'>#</a>
**`NEPSolver.newton`** &mdash; *Function*.



```
λ,v = newton([eltype],nep::NEP;[errmeasure,][tol,][maxit,][λ,][v,][c,][displaylevel,][armijo_factor=1,][armijo_max])
```

Applies Newton-Raphsons method on nonlinear equation with `(n+1)` unknowns. errmeasure is a function handle which specifies provides a procedure for error measure and termination. The iteration is continued until errmeausure is less than `tol`. `λ` and `v` are starting approximations. `c` is the orthogonalization vector.  If `c=0` the current approximation will be used for the orthogonalization. `armijo_factor` specifies if an Armijo rule should be applied, and its value specifies the scaling factor of the step length (per reduction step). The variable `armijo_max` specifies the maximum number of step length reductions.

**Example**

```julia-repl
julia> nep=nep_gallery("dep0");
julia> λ,v=newton(nep);
julia> minimum(svdvals(compute_Mder(nep,λ)))
1.6066157878930876e-16
```

**References**

  * Nichtlineare Behandlung von Eigenwertaufgaben, Z. Angew. Math. Mech. 30 (1950) 281-282.
  * A. Ruhe, Algorithms for the nonlinear eigenvalue problem, SIAM J. Numer. Anal. 10 (1973) 674-689

<a id='NEPSolver.augnewton' href='#NEPSolver.augnewton'>#</a>
**`NEPSolver.augnewton`** &mdash; *Function*.



```
Augmented Newton's method. Equivalent to newton() but works only with operations on vectors of length n, instead of n+1.
```

<a id='NEPSolver.resinv' href='#NEPSolver.resinv'>#</a>
**`NEPSolver.resinv`** &mdash; *Function*.



```
λ,v = resinv([eltype],nep::NEP;[errmeasure,][tol,][maxit,][λ,][v,][c,][displaylevel,][armijo_factor=1,][armijo_max,][linsolvecreator])
```

Applies residual inverse iteration method for nonlinear eigenvalue problems. linsolvecreator is a function which specifies how the linear system is created. See `newton()` for other parameters.

**Example**

```julia-repl
julia> nep=nep_gallery("qdep0");
julia> λ,v=resinv(Complex128,nep,λ=-2,v=ones(size(nep,1)))
julia> norm(compute_Mlincomb(nep,λ,v))
1.817030659827106e-14                   
```

**References**

  * A. Neumaier, Residual inverse iteration for the nonlinear eigenvalue problem, SIAM J. Numer. Anal. 22 (1985) 914-923

<a id='NEPSolver.quasinewton' href='#NEPSolver.quasinewton'>#</a>
**`NEPSolver.quasinewton`** &mdash; *Function*.



```
quasinewton([T=Complex128],nep,[errmeasure,][tol,][maxit,][λ,][v][ws][displaylevel][linsolvercreator,][armijo_factor,][armijo_max])
```

An implementation of the quasi-Newton approach referred to as quasi-Newton 2 in the reference. The method involves one linear system solve per iteration corresponding with the matrix $M(λ)$, where $λ$ is constant. The vector `ws` is a representation of the normalization, in the sense that $c^T=w_s^TM(λ)$, where all iterates satisfy $c^Tx_i=1$. See `newton()` for other parameters.

**Example**

```julia-repl
julia> nep=nep_gallery("pep0")
julia> λ,v=quasinewton(nep,v=ones(size(nep,1)));
julia> norm(compute_Mlincomb(nep,λ,v))/norm(v)
6.301479387102376e-15
```

**References**

  * Jarlebring, Koskela, Mele, Disguised and new Quasi-Newton methods for nonlinear eigenvalue problems, arxiv preprint: https://arxiv.org/abs/1702.08492

<a id='NEPSolver.mslp' href='#NEPSolver.mslp'>#</a>
**`NEPSolver.mslp`** &mdash; *Function*.



```
 mslp(nep,..)
```

Method of successive linear problems

<a id='NEPSolver.rfi' href='#NEPSolver.rfi'>#</a>
**`NEPSolver.rfi`** &mdash; *Function*.



```
rfi(nep,nept,[λ=0,][errmeasure=default_errmeasure,][tol=eps()*100,][maxit=100,][v=randn,][u=randn,][displaylevel=0,][linsolvecreator=default_linsolvecreator,])
```

Two-sided Rayleigh functional Iteration, as given as Algorithm 4 in  "Nonlinear Eigenvalue Problems: Newton-type Methods and Nonlinear Rayleigh Functionals", by Kathrin Schreiber.

<a id='NEPSolver.newtonqr' href='#NEPSolver.newtonqr'>#</a>
**`NEPSolver.newtonqr`** &mdash; *Function*.



```
Newton-QR method.
```

<a id='NEPSolver.implicitdet' href='#NEPSolver.implicitdet'>#</a>
**`NEPSolver.implicitdet`** &mdash; *Function*.



```
Implicit determinant method
```


<a id='Projection-methods-1'></a>

## Projection methods


```
NEPSolver.nlar
```


<a id='Arnoldi-type-methods-1'></a>

## Arnoldi type methods

<a id='NEPSolver.iar' href='#NEPSolver.iar'>#</a>
**`NEPSolver.iar`** &mdash; *Function*.



```
iar(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure=default_errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS])
```

**Infinite Arnoldi method**

Infinite Arnoldi method, as described in Algorithm 2 in  "A linear eigenvalue algorithm for the nonlinear eigenvalue problem", by Jarlebring, Elias and Michiels, Wim and Meerbergen, Karl.

<a id='NEPSolver.tiar' href='#NEPSolver.tiar'>#</a>
**`NEPSolver.tiar`** &mdash; *Function*.



```
tiar(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure=default_errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS])
```

**Tensor Infinite Arnoldi method**

Tensor Infinite Arnoldi method, as described in Algorithm 2 in  "The Waveguide Eigenvalue Problem and the Tensor Infinite Arnoldi Method", by Jarlebring, Elias and Mele, Giampaolo and Runborg, Olof.

<a id='NEPSolver.infbilanczos' href='#NEPSolver.infbilanczos'>#</a>
**`NEPSolver.infbilanczos`** &mdash; *Function*.



```
λv,V,U=infbilanczos([eltype],nep, nept,[linsolvecreator,][linsolvertcreator,][v,][u,][σ,][γ,][tol,][Neig,][errmeasure,][displaylevel,][maxit,][check_error_every])
```

Executes the Infinite Bi-Lanczos method on the problem defined by nep::NEP and nept::NEP. nep:NEP is the original nonlinear eigenvalue problem and nept::NEP is its (hermitian) transpose. v and u are starting vectors, σ is the shift and γ the scaling.  See `newton()` for other parameters.

**Example:**

```julia-repl
julia> nep=nep_gallery("dep0");
julia> A=get_Av(nep); fv=get_fv(nep);
julia> nept=SPMF_NEP([A[1]',A[2]',A[3]'],fv); # Create the transposed NEP
julia> λv,V=infbilanczos(nep,nept,Neig=3)
julia> norm(compute_Mlincomb(nep,λv[1],V[:,1]))
```

**References:**

  * The infinite bi-Lanczos method for nonlinear eigenvalue problems, S. W. Gaaf and E. Jarlebring, arxiv: 1607.03454, to appear in SIAM J. Scientific Computing


<a id='Gallery-1'></a>

# Gallery


WEP

