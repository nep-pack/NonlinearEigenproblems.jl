# Tutorial: Specifying linear solvers

Many of the NEP-solvers are based on solving linear systems of
the type
```math
M(λ)x=b.
```
In some methods the linear system matrices are the same, i.e., ``λ`` does not
change.
You can specify which numerical methods should be used to solve the linear system when you call a
NEP-solver. This tutorial illustrates this functionality,
and finally shows how you can specify your own method for linear systems.

## Built-in linear solvers

The linear solver is specified with the `linsolvercreator` keyword argument
in most NEP-solvers.
Let us contruct an example which we will solve with several methods.
The matrix ``M(λ)`` is sparse, and the nonlinearity is an exponential term:
```julia
using NonlinearEigenproblems, SparseArrays, LinearAlgebra;
n=100;
α=0.01;
A=spdiagm(0=>ones(n),1=>α*ones(n-1),-3=>α*ones(n-3));
B=spdiagm(0=>ones(n));
C=spdiagm(0=>(1:n)/n);
nep= SPMF_NEP([A,B,C],[s->one(s),s->s,s->exp(s)],align_sparsity_patterns=true);
λ0=-1.3; # Starting guess
```
Let us first solve it with the  [`resinv`](@ref) method, using the default solver for the linear system:
```julia-repl
julia> (λ,x)=resinv(nep,λ=λ0,v=ones(n),logger=1,tol=1e-16);
Precomputing linsolver
iter 1 err:0.0066371687626530325 λ=-1.3 + 0.0im
iter 2 err:0.005517924619612248 λ=-1.175478914232863 + 0.0im
iter 3 err:0.002478390282482615 λ=-1.237914206446667 + 0.0im
iter 4 err:0.0007175125397164684 λ=-1.2715354474842264 + 0.0im
...
iter 68 err:1.9815836266850424e-16 λ=-1.2845622481786096 + 0.0im
iter 69 err:1.2398848173585647e-16 λ=-1.28456224817861 + 0.0im
iter 70 err:8.341032972349128e-17 λ=-1.2845622481786103 + 0.0im

```
We will carry out some timing experiments, so let's
use the [`BenchmarkTools`](https://github.com/JuliaCI/BenchmarkTools.jl)-package and switch off printouts in
the NEP-solver:
```julia-repl
julia> using BenchmarkTools
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),tol=1e-16);
  8.170 ms (32457 allocations: 10.99 MiB)
```
The linear system that has to be solved in every iteration
in `resinv` has a constant system matrix, and
therefore a prefactorization (typically an LU-factorization)
is useful. This is done with the [`FactorizeLinSolverCreator`](@ref),
which is actually the default behaviour, so we
get no substantial difference when we specify
a creator if the type `FactorizeLinSolverCreator`.
```julia-repl
julia> creator=FactorizeLinSolverCreator();
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator,tol=1e-16);
  8.104 ms (32447 allocations: 10.98 MiB)
```
If we do not want to use a prefactorization, you can specify
[`BackslashLinSolverCreator`](@ref) as your creator object.
```julia-repl
julia> creator=BackslashLinSolverCreator();
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator,tol=1e-16);
  19.985 ms (36932 allocations: 22.31 MiB)
```
This does not use a prefactorization and is therefore slower.


The above approach corresponded to direct methods for linear systems.
You can also use iterative methods, e.g., the GMRES-method.
The GMRES-method is available in the [`GMRESLinSolverCreator`](@ref)
function. Iterative methods in general need [preconditioners](https://en.wikipedia.org/wiki/Preconditioner).
We continue the example and use a diagonal preconditioner:
```julia-repl
julia> D0=(Diagonal(compute_Mder(nep,λ0))); # Preconditioner
julia> creator=GMRESLinSolverCreator(Pl=D0, tol=1e-2);
```
All the keyword arguments in the call `GMRESLinSolverCreator`
are passed to
[`gmres!`](https://juliamath.github.io/IterativeSolvers.jl/stable/linear_systems/gmres/).
Hence, the `tol` here  specifies a termination criteria for the GMRES-method,
and `Pl` specifies the left preconditioner, in this case just a diagonal matrix.
```julia-repl
julia> (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator,logger=1,tol=1e-16);
Precomputing linsolver
iter 1 err:0.0066371687626530325 λ=-1.3 + 0.0im
iter 2 err:0.005516840431746622 λ=-1.175478914232863 + 0.0im
iter 3 err:0.002478456565315529 λ=-1.2379013065364441 + 0.0im
iter 4 err:0.0007137260829914206 λ=-1.271604846382212 + 0.0im
...
iter 69 err:1.4948517196633335e-16 λ=-1.2845622481786099 + 0.0im
iter 70 err:1.0136830543996086e-16 λ=-1.28456224817861 + 0.0im
iter 71 err:6.468407765141646e-17 λ=-1.2845622481786103 + 0.0im

```
The printout reveals that we need one more more iteration, than with a direct method.
In terms of computation time, this approach can however still be competitive:
```julia-repl
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator,tol=1e-16);
  10.912 ms (49183 allocations: 18.95 MiB)
```

## Your own linear solver

There are many ways to solve linear systems in Julia, e.g.,
by using package such as [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl),
[Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl)
or [KrylovMethods.jl](https://github.com/JuliaInv/KrylovMethods.jl).
These are not natively supported by NEP-PACK,
but due to the extendability of the `LinSolverCreator`-objects
specified above, you can still use them.
We illustrate the extendability by creating a linear solver
based on solving a [Schur complement](https://en.wikipedia.org/wiki/Schur_complement).
The following helper-function
for the Schur complement solve will be used later.
```julia
function schur_complement_lin_solve(AA,b,n0)
  A=AA[1:n0,1:n0];
  B=AA[1:n0,(n0+1):end];
  C=AA[(n0+1):end,1:n0];
  D=AA[(n0+1):end,(n0+1):end];
  S=D-C*(A\B); # Schur complement
  b1=b[1:n0]; b2=b[(n0+1):end];
  Ainvb1=A\b1; Sinvb2=S\b2;
  # Formula for the linear solve:
  x1=A\(b1+(B*(S\(C*(Ainvb1)))))-A\(B*(Sinvb2))
  x2=-S\(C*(Ainvb1))+Sinvb2;
  return [x1;x2];
end
```

Julia's efficiency stems partially from the extensive use of types.
We need to define new types to specify our own linear solver
and integrate it with NEP-PACK.

```julia
struct MyLinSolverCreator <: LinSolverCreator; end
struct MyLinSolver <: LinSolver;
  mynep
  myλ
end
```
NEP-solvers call the function [`create_linsolver`](@ref)(creator,nep,λ),
which should return a linear solver. We need to overload this function
for our own creator-type.
In general, this is to allow precomputation.
However, in this example we do not have any precomputations and
thus just return an instance of `MyLinSolver`.
```julia
import NonlinearEigenproblems.create_linsolver # Needed since we want overload it
function create_linsolver(::MyLinSolverCreator,nep,λ)
   return MyLinSolver(nep,λ);
end
```
The rest of the implementation of the solver goes into the function `lin_solve`, where we
utilize our function `schur_complement_lin_solve` from above.
```julia
import NonlinearEigenproblems.LinSolvers.lin_solve # Needed since we want overload it
function lin_solve(solver::MyLinSolver,b::Vector;tol=eps())
   n0=10;
   return schur_complement_lin_solve(compute_Mder(solver.mynep,solver.myλ),b,n0)
end
```
You can now solve the problem by passing a creator object `MyLinSolverCreator()` to a
NEP-solver, e.g., [`augnewton`](@ref):
```julia-repl
julia> dep=nep_gallery("dep0",50);
julia> creator=MyLinSolverCreator();
julia> augnewton(dep,λ=1,v=ones(size(dep,1)),logger=1,linsolvercreator=creator);
iter 1 err:0.10615052208736536 λ=1.0 + 0.0im
iter 2 err:0.04682362994161844 λ=3.004830719411172 + 0.0im
iter 3 err:0.08148964717804698 λ=0.213140384062811 + 0.0im
iter 4 err:0.03955381142282053 λ=0.47667949368248896 + 0.0im
iter 5 err:0.06584371583464586 λ=2.985447356041631 + 0.0im
iter 6 err:0.02262384918568079 λ=3.722422057973499 + 0.0im
iter 7 err:0.0036373167693678717 λ=3.389502913018821 + 0.0im
iter 8 err:0.00704620404184537 λ=3.2745554693864496 + 0.0im
iter 9 err:0.0009450496517445758 λ=3.1652287152386758 + 0.0im
iter 10 err:2.138372017573122e-5 λ=3.187725547526568 + 0.0im
iter 11 err:1.0100960591678548e-8 λ=3.188230946035159 + 0.0im
iter 12 err:2.801564990446382e-15 λ=3.1882313460682705 + 0.0im
```

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_LINSOLVE)
