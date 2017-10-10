
# NEPCore Documentation



```@docs
NEP
```

```@docs
SPMF_NEP
```

```@docs
SPMF_NEP(AA,fii,Schur_fact=false)
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

# Compiling the documentation

Compile this documentation site by running:
```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ julia --color=yes make.jl
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ cat build/index.md |  perl -pe 's/([^\$])\$([^\$]+)\$/$1\\\\($2\\\\)/g'  > build/index2.md
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ mkdocs build --clean
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ firefox site/index2/index.html
```
The perl-expression is a dirty trick to make inline math display properly. 

More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)

