# Tutorial: Solving NEP defined in MATLAB

## A problem defined in MATLAB

MATLAB is the de-facto standard language for many tasks in scientific
computing.
If you have a NEP defined in MATLAB, you can quite easily
use the NEP-solvers of this package. Below is a description
of two ways to solve nonlinear eigenvalue
problems defined in MATLAB.
There is a cost in terms of efficiency to define your problem in MATLAB,
due to overhead associated with communication between the
MATLAB and Julia processes. Very large scale problems are recommended to be defined
directly in Julia.

!!! note
    To work with NEPs defined in MATLAB you need to have MATLAB installed on
    your computer. We use the [MATLAB interoperability package](https://github.com/JuliaInterop/MATLAB.jl)
    to link Julia execution with MATLAB execution.

Suppose you have the following NEP in MATLAB
```math
M(\lambda)=A_0+\lambda A_1+\exp(\lambda A_2),
```
where ``A_1,A_2,A_3`` are martices and ``\exp``
the matrix exponential.
The problem can be defined in MATLAB as follows.
This is the contents of the file `compute_derivative_k.m`:
```matlab
function Z=compute_derivative_k(s,k)
     randn('seed',0);
     n=10;
     A0=randn(n,n); A1=randn(n,n);
     Z=zeros(n,n);
     if (k==0)
         Z=A0+s*A1;
     end
     if (k==1)
         Z=A1;
     end
     Z=Z+(A1^k)*expm(s*A1);
end
```
The function which computes derivative `k` evaluted in the point `s`.
We assume in the following that the file `compute_derivative_k.m`
is located in the current directory.

## Approach 1: Implementation in NEP-PACK (using `Mder_NEP`)

The easiest way to create a NEP which is only defined by
its derivative computation is by the helper type
[`Mder_NEP`](@ref).
```julia
using NonlinearEigenproblems, MATLAB
nep=Mder_NEP(10,(s,der) -> mat"compute_derivative_k($s,double($der))");
```
The NEP can now be approached with many of the methods in the
package, e.g., with a contour integral method:
```julia-repl
julia> (λ,V)=contour_beyn(nep, radius=0.6, k=8);
julia> λ
2-element Array{Complex{Float64},1}:
 0.1711954796771912 - 6.401495587242332e-15im
 0.1547216302712358 - 0.16631220583083045im   
```
The first argument of the `Mder_NEP` instantiation
is the size of the NEP.
The instantiation of the `Mder_NEP` creates a NEP-object
only defined by its matrix derivative functions,
given in the call-back function specified by
the second argument.
In this case, the the function calls a
MATLAB process (running in the background completely hidden
from the Julia user) and requests a execution of
`compute_derivate_k` with the given arguments. After executing the
MATLAB-call, the MATLAB-process sends the matrix back to Julia.
In other words, we have coupled the derivative
computation of the NEP with a call to MATLAB.
More precisely, every call to the [`compute_Mder`](@ref) function
leads to a call to the created MATLAB function. Compare:
```julia-repl
julia> compute_Mder(nep,0.1+0.2im)
10×10 Array{Complex{Float64},2}:
   2.15543-0.101289im     -1.4933-0.0817057im  …    0.220658-0.313894im    0.371533+0.544535im
  0.751339+0.213974im    0.532557+0.359256im         0.58971-0.135805im    -2.65352-0.308557im
 -0.177809-0.383021im     1.46944+0.374974im        0.965219+0.140168im    0.979515-0.603281im
  0.204312-0.300014im   -0.576669+0.630099im         1.94032+0.43922im    -0.803364+0.652243im
 -0.378807+0.511258im    0.590185+0.0812181im      -0.381726-0.138047im     1.30005+0.562022im
   1.58041-0.266624im   -0.347895+0.292268im   …  -0.0826327+0.155039im   -0.691359+0.340299im
  0.105255+0.0940046im  -0.338352-0.443379im        0.546516-0.062307im    0.445814-0.648929im
   1.47467-0.646341im    -1.36607-0.195403im         1.02226-0.0228401im    1.14153+0.576545im
  0.182295-0.143594im    -1.06099+0.492347im        -0.41297-0.409332im   -0.322912-0.219094im
  0.658504-0.190844im     1.21896+0.280606im       -0.563413-0.073228im     1.25092-0.418521im
```
with the MATLAB call:
```matlab
>> M=compute_derivative_k(0.1+0.2i,0);
>> M(1:3,1:3) % I don't want to see the whole matrix
ans =
   2.1554 - 0.1013i  -1.4933 - 0.0817i   0.1131 + 0.0836i
   0.7513 + 0.2140i   0.5326 + 0.3593i  -1.0211 - 0.5134i
  -0.1778 - 0.3830i   1.4694 + 0.3750i   0.2180 + 0.1720i
```
You can verify that the output of the
call to the [`contour_beyn`](@ref)-method is a solution
directly in MATLAB:
```matlab
>> s = 0.1547216302712358 - 0.16631220583083045i; % copied from the output above (remember: 1im -> 1i)
>> M=compute_derivative_k(s,0);
>> min(svd(M)) % Matrix is singular if s is a solution
ans =
   1.5239e-15
```

NEP-objects in NEP-PACK are defined from compute-functions (as
we describe in [NEPTypes](types.md)) and in this case we only defined
the derivative computation function `compute_Mder`. Note that
the `Mder_NEP`-type
provides default implementations of
[`compute_Mlincomb`](@ref) as well as [`compute_MM`](@ref) (by wrapping
calls to `compute_Mder`) in a
way that is hidden from the
user, such that we can still use algorithms
based on those compute functions.
More efficiency can be obtained if these compute
functions are also implemented, e.g., by
a different MATLAB-function.

## Approach 2: Implementation in NEP-PACK (using new type)

We illustrate the extendability of the package
by defining our own type, which again
uses the MATLAB-package in the background.

!!! note
    The process is also described in the
    [BEM tutorial](bemtutorial.md#Implementation-in-NEP-PACK-using-the-Mder_NEP-type-1).

The size is hardcoded in this example, so we can
define a new type of the specific size:
```julia
struct MATLABNEP <: NEP
end
Base.size(nep::MATLABNEP) = (10,10)
Base.size(nep::MATLABNEP,::Int) = 10
```
Initiate the MATLAB package and prepare to integrate with NEP-PACK:
```julia
using MATLAB; # requires MATLAB to be installed
mat"addpath('.')" # Add path to your m-file
import NonlinearEigenproblems.compute_Mder;
import NonlinearEigenproblems.compute_Mlincomb;
import NonlinearEigenproblems.compute_Mlincomb_from_Mder;
```

In this example, the problem is only provided by a function
to compute derivatives of `M`, which
we specify by defining a  `compute_Mder` function.
We also specify that linear combinations of derivatives should
be computed by calling `compute_Mder` in the naive way:
```julia
function compute_Mder(::MATLABNEP,s::Number,der::Integer=0)
    return mat"compute_derivative_k(double($s),double($der))"
end
compute_Mlincomb(nep::MATLABNEP,λ::Number,V::AbstractVecOrMat, a::Vector) = compute_Mlincomb_from_Mder(nep,λ,V,a)
compute_Mlincomb(nep::MATLABNEP,λ::Number,V::AbstractVecOrMat) = compute_Mlincomb(nep,λ,V, ones(eltype(V),size(V,2)))
```
Now you can instantiate the NEP and use your favorite NEP-solver,
in this case we use [`newtonqr`](methods.md#NonlinearEigenproblems.NEPSolver.newtonqr).
```julia-repl
julia> nep=MATLABNEP();
julia> (λ,v)=newtonqr(nep,λ=-3,logger=1,maxit=30,v=ones(10));
iter 1 err:1.033593309412195 λ=-3.0 + 0.0im
iter 2 err:0.3059246224011592 λ=0.83641207310996 + 0.0im
iter 3 err:0.6000405834026614 λ=-1.7728647881500432 + 0.0im
iter 4 err:0.07375061614602237 λ=-0.7800560594951582 + 0.0im
iter 5 err:0.0093516562758152 λ=-0.8707521093182906 + 0.0im
iter 6 err:8.954564848847882e-5 λ=-0.8840785307305598 + 0.0im
iter 7 err:7.446596609664408e-9 λ=-0.88420751056806 + 0.0im
iter 8 err:1.0942739518352825e-15 λ=-0.884207521294992 + 0.0im
```
The residual is small and we have a solution
```julia-repl
julia> using LinearAlgebra
julia> norm(compute_Mlincomb(nep,λ,v))/norm(v)
1.0942739518352825e-15
```


![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_MATLAB1)
