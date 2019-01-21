# Tutorial: Solving NEP defined in MATLAB

## A problem defined in MATLAB

MATLAB is a de-facto standard for many tasks in scientific
computing.
If you have a NEP defined in MATLAB, you can quite easily
use the NEP-solvers of this package. Below is a description
of one way to interface with MATLAB. The example
illustrates the principle at the cost of some efficiency.

Suppose you have the following NEP in MATLAB
```math
M(\lambda)=A_0+\lambda A_1+\exp(\lambda A_2).
```
The problem can be defined in MATLAB as follows.
This is the contents of the file `compute_derivative_k.m`

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

## Implementation in NEP-PACK

We define a new type representing our MATLAB-NEP.
The size is hardcoded in this example.
```julia
struct MATLABNEP <: NEP
end
Base.size(nep::MATLABNEP) = (10,10)
Base.size(nep::MATLABNEP,::Int) = 10
```
Initiate the MATLAB package and prepare to integrate with NEP-PACK:
```julia
julia> using MATLAB; # requires MATLAB to be installed
julia> mat"addpath('.')" # Add path to your m-file
julia> import NonlinearEigenproblems.compute_Mder;
julia> import NonlinearEigenproblems.compute_Mlincomb;
julia> import NonlinearEigenproblems.compute_Mlincomb_from_Mder;
```
NEP-objects in NEP-PACK are defined from compute-functions (as
we describe in [NEPTypes](types.md)) and we need to define
the derivative computation function, which calls the MATLAB-code.
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
```julia
julia> nep=MATLABNEP();
julia> (λ,v)=newtonqr(nep,λ=-3,displaylevel=1,maxit=30,v=ones(10))
Iteration: 1 errmeasure: 1.0335933094121779
Iteration: 2 errmeasure: 0.305924622401145
Iteration: 3 errmeasure: 0.6000405833925101
Iteration: 4 errmeasure: 0.07375061613894424
Iteration: 5 errmeasure: 0.009351656273646538
Iteration: 6 errmeasure: 8.954564844507815e-5
Iteration: 7 errmeasure: 7.446596374256243e-9
Iteration: 8 errmeasure: 1.8095439571351245e-15
(-0.8842075212949918 + 0.0im, Complex{Float64}[0.544936+0.0im, 0.641218+0.0im, 0.089366+0.0im, -0.0975048+0.0im, 0.133397+0.0im, 1.0+0.0im, -0.836009+0.0im, -0.00753176+0.0im, 0.270149+0.0im, -0.664448+0.0im], [0.354722, -0.0659026, -0.465767, 0.079273, -0.524316, -0.372411, -0.0129146, -0.386585, -0.140157, 0.252488])
```
The residual is small and we have a solution
```julia
julia> norm(compute_Mlincomb(nep,λ,v))
3.111596859559977e-15
```
![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_MATLAB1)
