
<a id='NEPCore-Documentation-1'></a>

# NEPCore Documentation

<a id='NEPCore.NEP' href='#NEPCore.NEP'>#</a>
**`NEPCore.NEP`** &mdash; *Type*.



```
abstract NEP
```

NEP represents a nonlinear eigenvalue problem

<a id='NEPTypes.SPMF_NEP' href='#NEPTypes.SPMF_NEP'>#</a>
**`NEPTypes.SPMF_NEP`** &mdash; *Type*.



```
type SPMF_NEP <: AbstractSPMF
```

An SPMF_NEP is defined by the sum the sum $Σ_i A_i f_i(λ)$, where i = 0,1,2,..., all of the matrices are of size n times n and f_i is a function. In particular, it must be possible to evaluate f_i with a matrix argument
Constructor: SPMF_NEP(AA,fii,Schur_fact = false) where AA is an array of the matrices A_i and  fii is an array of the funtion f_i. Set $Schur_fact = true$ if you want to pre-factorize the matrices in the call of $compute_MM(...)$.

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


<a id='Compiling-the-documentation-1'></a>

# Compiling the documentation


Compile this documentation site by running:


```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ julia --color=yes make.jl
```

