
# NEPCore Documentation



```@docs
NEP
```

```@docs
SPMF_NEP
```

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

# Compiling the documentation

Compile this documentation page by running:
```
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ julia --color=yes make.jl
jarl@bjork:~/jobb/src/nep-pack-alpha/docs$ mkdocs build --clean
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
jarl@bjork:/tmp/nep-pack-alpha$ cp -r $DOCSDIR/* .
jarl@bjork:/tmp/nep-pack-alpha$ git add *
jarl@bjork:/tmp/nep-pack-alpha$ git commit .
jarl@bjork:/tmp/nep-pack-alpha$ git push
```


More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)




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

WEP


