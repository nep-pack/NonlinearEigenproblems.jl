# Tutorial: Using NEP-PACK from python

## PyJulia

The previous tutorial illustrated how a NEP defined
in python code can be solved with NEP-PACK.
Python is currently a more mature language than Julia,
and there are considerable packages and features
in python not present in Julia. If you need these features,
it may be more convenient to call NEP-PACK
from  python, rather than calling python code from julia.

The python package [PyJulia](https://github.com/JuliaPy/pyjulia)
gives us that possibility. The installation of PyJulia on Ubuntu linux is simple:
```
$ python3 -m pip install julia # Only necessary first time you run it
...
$ python3
>>> from julia.api import Julia
>>> jl = Julia(compiled_modules=False) # compilation flag necessary on ubuntu
>>> julia.install()               # Only necessary first time you run it
>>> from julia import Base
>>> Base.MathConstants.golden  # Julia's definition of golden ratio
1.618033988749895
```


## Using PyJulia and NEP-PACK

The [`Mder_NEP`](@ref)-function provides a convenient
way to define NEPs by only
using a function that computes the matrix ``M(λ)``
and its derivatives.
Let us first define a function which does that in python. We consider
the problem
```math
M(λ)=\begin{bmatrix}3&2\newline3&-1\end{bmatrix}+
λ\begin{bmatrix}0&2\newline0&1\end{bmatrix}+
e^{0.5 λ}\begin{bmatrix}1&1\newline1&1\end{bmatrix}
```
and implement it with this python code:
```python
import numpy as np;
import cmath as m;
def my_compute_M(s,der):
    """Compute the matrix M^{(k)}(s) for a given eigenvalue approximation and derivative k"""
    A=np.matrix('3 2; 3 -1');  B=np.matrix('0 2; 0 1');   C=np.matrix('1 1; 1 1');
    tau=0.5;
    M=pow(tau,der)*m.exp(tau*s)*C
    if (der==0):
        M=M+A+s*B;
    elif (der==1):
        M=M+B;
    return M
```
An evaluation of the matrix function can be done by the call:
```
>>> my_compute_M(0.3,0)
matrix([[4.16183424+0.j, 3.76183424+0.j],
        [4.16183424+0.j, 0.46183424+0.j]])
```
We instantiate a new NEP based with `Mder_NEP` which first must be imported
```
>>> from julia.NonlinearEigenproblems import Mder_NEP
>>> n=2; # Size of the problem
>>> nep=Mder_NEP(2,my_compute_M);
```
and we can apply most of our solvers to this problem by first importing the corresponding function, in this case we use [`contour_beyn`](@ref).
```
>>> from julia.NonlinearEigenproblems import contour_beyn;
>>> sol=contour_beyn(nep,logger=1,neigs=1,radius=3)
Computing integrals
NonlinearEigenproblems.NEPSolver.MatrixTrapezoidal: computing G...
NonlinearEigenproblems.NEPSolver.MatrixTrapezoidal: summing terms........................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
Computing SVD prepare for eigenvalue extraction  p=1
Computing eigenvalues
Computing eigenvectors
>>>
```
We can verify that we computed a solution as follows
```
>>> s=sol[0][0]; v=sol[1]
>>> my_compute_M(s,0)*v
matrix([[1.71634841e-17-1.59872116e-14j],
        [9.55210099e-17-3.99680289e-15j]])
>>> from numpy.linalg import norm
>>> norm(my_compute_M(s,0)*v)
1.6479526251408437e-14
```
Note that in order to obtain better efficiency for
large-scale problems, and reduce overhead,
you may want to use [`Mder_Mlincomb_NEP`](@ref),
as described in the [previous tutorial](tutorial_call_python.md).


![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_PYTHON2)
