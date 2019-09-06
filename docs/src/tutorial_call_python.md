# Tutorial: Solving NEP defined in Python

## A problem defined in Python

Julia is a great programming language,
but your problem may not be easy to define in Julia code, e.g., for legacy reasons.
Don't let that prevent you from using the package.
We now show how a problem defined in Python can be solved
with NEP-PACK.

One of the advantages of the Julia language is that it
is reasonably easy to interface with code written in
other languages. In this tutorial we work with Python, and the two
following tutorials we interface [MATLAB](tutorial_matlab1.md) and
[fortran](tutorial_fortran1.md).

!!! note
    To work with NEPs defined in [Python](https://www.python.org/) you
    need to have Python installed on
    your computer. With the pakcage [PyCall](https://github.com/JuliaPy/PyCall.jl)
    it is possible to let Julia control execution of Python code.

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

## Interfacing Python code

We first initiate `PyCall` as follows. Note that the `pushfirst!` command
is needed, otherwise the module defined in the file `mynep.py` we gave
above will not be found.
(`PyCall` does not include the current directory in the module search path by default.)

```julia
using PyCall;
pushfirst!(PyVector(pyimport("sys")."path"), "");
mynep = pyimport("mynep")
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

## Implementation in NEP-PACK (using `Mder_Mlincomb_NEP`)
We can now use the Python interface to define a NEP in Julia.
The type [`Mder_Mlincomb_NEP`](@ref) is a special type made for this situation.
The required inputs are the size, called `n`; a function to compute
$M(λ)$, called `fder`; and a function
to compute $\sum_{i=1}^kM^{(k)}(λ)x_i$, called `flincomb`.
The extra `0` passed in the definition defines that $M(λ)$ is available,
but no higher derivatives.
```julia
using NonlinearEigenproblems
n=2;
fder = (λ,der) -> mynep.compute_M(complex(λ));
flincomb =  (λ,X) -> mynep.compute_Mlincomb(complex(λ),complex(reshape(X,size(X,1),size(X,2))));
nep=Mder_Mlincomb_NEP(n,fder,0,flincomb);
```
We can compare the Python call with the NEP-PACK call
```julia-repl
julia> compute_Mder(nep,3+3im)
2×2 Array{Complex{Float64},2}:
 -18.8845+2.83447im  -17.8845+2.83447im
 -16.8845+2.83447im  -12.8845+5.83447im
```
We continue by computing some eigenvalues of the the NEP using the
Infinite Arnoldi method ([`iar`](@ref)).
```julia-repl
julia> (λ,v)=iar(nep,v=[1;1],σ=1,logger=0,neigs=3);
julia> λ
3-element Array{Complex{Float64},1}:
  0.6748316143423988 + 7.336803319821954e-19im
 0.11742590291190791 - 3.649946317867008im    
 0.11742590291191168 + 3.6499463178670144im  
```
We can verify that we actually computed solutions as follows:
```julia-repl
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
1.106424240899132e-15
```

## Implementation in NEP-PACK  (using new type)
The previous implementation utilizes the convenience type `Mder_Mlincomb_NEP`,
and solves the problem in a satisfactory way. Nevertheless, to illustrate more
of the inner workings of NEP-PACK we solve the problem in a second way,
by defining our own type.
The first thing we need to do is to define the `size`-function, which is
hard-coded in this example.
```julia
import NonlinearEigenproblems.size # We will overload these functions
import NonlinearEigenproblems.compute_Mlincomb;
import NonlinearEigenproblems.compute_Mder;
struct PyNEP <: NEP # Set up a dummy type for our specific NEP
end
size(::PyNEP) = (2,2) # Trivial function definitions
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
julia> pynep=PyNEP();
julia> compute_Mder(pynep,3+3im)
2×2 Array{Complex{Float64},2}:
 -18.8845+2.83447im  -17.8845+2.83447im
 -16.8845+2.83447im  -12.8845+5.83447im
```
The behavior is the same as above.
Since a NEP-object is defined by its compute functions,
we can now use many NEP-solvers to solve this problem.
We again use [`iar`](@ref):
```julia-repl
julia> (λ2,v2)=iar(pynep,v=[1;1],σ=1,logger=0,neigs=3);
julia> λ2
3-element Array{Complex{Float64},1}:
  0.6748316143423988 + 7.336803319821954e-19im
 0.11742590291190791 - 3.649946317867008im    
 0.11742590291191168 + 3.6499463178670144im   
```
We can compare with the eigenvalues computed above and, again,
verify that we actually computed solutions as follows:
```julia-repl
julia> norm(compute_Mlincomb(pynep,λ2[1],v2[:,1]))
1.106424240899132e-15
```
Residual is almost zero, so we have a solution.

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_PYTHON1)
