# Tutorial: Solving NEP defined in Python

## A problem defined in Python

Julia is a great programming language,
but your problem may not be easy to define in Julia code, e.g., for legacy reasons.
Don't let that prevent you from using the package.
We now show how a problem defined in python can be solved
with NEP-PACK.

The following python code correspond to the NEP
```math
M(λ)=\begin{bmatrix}1&2\newline3&4\end{bmatrix}+
λ\begin{bmatrix}0&0\newline0&1\end{bmatrix}+
e^λ\begin{bmatrix}1&1\newline1&1\end{bmatrix}
```
The code has two functions:
one that can compute an evaluation of $M(λ)$ and
one that can form a linear combination of derivatives
```math
  \sum_{i=1}^kM^{(k)}(λ)x_i.
```
Put a file  `mynep.py`  in your current directory with the following contents:
```python
import numpy as np;
import cmath as m;
def compute_M(s):
    """Compute the matrix M(s) for a given eigenvalue approximation"""
    A=np.matrix('1 2; 3 4');  B=np.matrix('0 0; 0 1');   C=np.matrix('1 1; 1 1');
    M=A+s*B+m.exp(s)*C
    return M

def compute_Mlincomb(s,X):
    """Compute the linear combination of derivatives for value s"""
    A=np.matrix('1 2; 3 4');  B=np.matrix('0 0; 0 1');   C=np.matrix('1 1; 1 1');

    X=np.matrix(X) # Explicitly convert to matrix
    z=np.zeros((2,1));
    # Zeroth derivative
    z=z+A*X[:,0]
    z=z+B*(s*X[:,0])
    z=z+C*(m.exp(s)*X[:,0])

    # First derivative
    if (np.size(X,1)>1):
        z=z+B*(X[:,1])+C*(m.exp(s)*X[:,1])
    # Higher derivatives
    if (np.size(X,1)>1):
        for k in range(2,np.size(X,1)):
            z=z+C*(m.exp(s)*X[:,k])
    return z
```

## Implementation in NEP-PACK

One of the advantages of the Julia language is that it
is reasonably easy to interface with code written in
other langauges. The Julia package [PyCall](https://github.com/JuliaPy/PyCall.jl)
simplifies the use of Python code and Julia code.


We first initiate `PyCall` as follows. Note that the `pushfirst!` command
is needed, otherwise the module defined in the file `mynep.py` we gave
above will not be found. (`PyCall` does not include the current directory in the module search path by default.)

```julia
using PyCall
pushfirst!(PyVector(pyimport("sys")["path"]), "")
local mynep
@pyimport mynep as mynep
```
This gives us direct access to the `compute_M`
and `compute_Mlincomb` functions from python, e.g.,
if we want to evaluate $M(3+3i)$ we run this code
```julia-repl
julia> mynep.compute_M(3+3im)
2×2 Array{Complex{Float64},2}:
 -18.8845+2.83447im  -17.8845+2.83447im
 -16.8845+2.83447im  -12.8845+5.83447im
```

We now just need to define a NEP which calls these routines.
It is achieved by defining a new NEP-PACK type, for
which we have define the `size`-function, which is
hard-coded in this example.

```julia
using NonlinearEigenproblems
import NonlinearEigenproblems.size
import NonlinearEigenproblems.compute_Mlincomb;
import NonlinearEigenproblems.compute_Mder;
struct PyNEP <: NEP # Set up a dummy type for our specific NEP
end
size(::PyNEP) = (2,2)
size(::PyNEP,::Int) = 2
```
As explained in [NEPTypes](types.md), a NEP is defined by
its compute functions. Here is how you define two compute functions
that call our python-defined NEP:
```julia
function compute_Mder(::PyNEP,s::Number,der::Integer=0)
    if (der>0)
        error("Higher derivatives not implemented");
    end
    return mynep.compute_M(complex(s)); # Call python
end
function compute_Mlincomb(::PyNEP,s::Number,X::AbstractVecOrMat)
    XX=complex(reshape(X,size(X,1),size(X,2))) # Turn into a matrix
    return mynep.compute_Mlincomb(complex(s),XX); # Call python
end
```
We now create an object of our newly created type and we can access the
NEP with the NEP-PACK specific compute functions:
```julia-repl
julia> pnep=PyNEP();
julia> compute_Mder(pnep,3+3im)
2×2 Array{Complex{Float64},2}:
 -18.8845+2.83447im  -17.8845+2.83447im
 -16.8845+2.83447im  -12.8845+5.83447im
```

## Solving the NEP

Since a NEP-object is defined by its compute functions,
we can now use many NEP-solvers to solve this problem.
Here we use IAR:
```julia-repl

julia> (λv,vv)=iar(pnep,v=[1;1],σ=1,logger=1,neigs=3);
Iteration:1 conveig:0
Iteration:2 conveig:0
Iteration:3 conveig:0
Iteration:4 conveig:0
Iteration:5 conveig:0
Iteration:6 conveig:0
Iteration:7 conveig:0
....
Iteration:26 conveig:1
Iteration:27 conveig:1
Iteration:28 conveig:1
julia>
```
We can verify that we actually computed solutions as follows:
```julia-repl
julia> λ=λv[1]; # Take the first eigenpair
julia> v=vv[:,1]
2-element Array{Complex{Float64},1}:
 -0.7606536306084172 + 4.723354443026557e-18im
   0.568748796395112 + 1.8318449036023953e-19im
julia> A=[1 2 ; 3 4];
julia> B=[0 0 ; 0 1];
julia> C=[1 1 ; 1 1];
julia> r=A*v+λ*B*v+exp(λ)*C*v;
2-element Array{Complex{Float64},1}:
 -3.3306690738754696e-16 + 1.4448222154182884e-17im
 -1.0547118733938987e-15 + 2.4802198512062408e-17im
```
Residual is almost zero, so we have a solution.

Note: The above functionality can also be achieved with  `Mder_NEP` in the development version of NEP-PACK

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_PYTHON1)
