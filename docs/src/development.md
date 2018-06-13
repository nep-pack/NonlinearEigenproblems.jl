# Developer info


# Compiling the documentation

Compile this documentation page by running:
```
jarl@bjork:~/src/NonlinearEigenproblems.jl/docs$ julia --color=yes make.jl &&  mkdocs build --clean
jarl@bjork:~/src/NonlinearEigenproblems.jl/docs$ firefox site/index.html
```
If you want this to appear on our documentation page
[https://nep-pack.github.io/NonlinearEigenproblems.jl/](https://nep-pack.github.io/NonlinearEigenproblems.jl/)
you need to push it to the `gh-branch`, e.g.,  by running
```
jarl@bjork:~/src/NonlinearEigenproblems.jl/docs$ export DOCSDIR=`pwd`
jarl@bjork:~/src/NonlinearEigenproblems.jl/docs$ cd /tmp
jarl@bjork:/tmp$ git clone -b "gh-pages" git@github.com:nep-pack/NonlinearEigenproblems.jl.git
jarl@bjork:/tmp$ cd NonlinearEigenproblems.jl
jarl@bjork:/tmp/NonlinearEigenproblems.jl$ cp -r $DOCSDIR/site/* .
jarl@bjork:/tmp/NonlinearEigenproblems.jl$ git add *;  git commit . -m "refresh docs"; git push
```


More information about `Documenter.jl`: [here](https://juliadocs.github.io/Documenter.jl/v0.1.3/man/guide/#Package-Guide-1)

