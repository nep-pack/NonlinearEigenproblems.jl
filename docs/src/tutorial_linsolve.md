# Tutorial: Your own linear solver

Many of the NEP-solvers are based on solving linear systems of
the type
```math
M(λ)x=b.
```
In some methods the linear system matrices are the same, i.e., `λ` does not
change.
You can specify which numerical methods should be used when you call a
NEP-solver. This tutorial illustrates this functionality,
and finally shows how you can specify your own method for linear systems.

## Built-in linear solvers

The linear solver is specified with the `linsolvercreator` keyword argument
in most NEP-solvers.
Let us contruct an example which we will solve with several methods.
It is a sparse matrix with an exponential term:
```julia
julia> using SparseArrays,LinearAlgebra;
julia> n=100;
julia> α=0.01;
julia> A=spdiagm(0=>ones(n),1=>α*ones(n-1),-3=>α*ones(n-3));
julia> B=spdiagm(0=>ones(n));
julia> C=spdiagm(0=>(1:n)/n);
julia> nep= SPMF_NEP([A,B,C],[s->one(s),s->s,s->exp(s)],align_sparsity_patterns=true);
```
Let us first solve it with the  [`resinv`](@ref) method.
```julia
julia> λ0=-1.2; # Starting guess
julia> (λ,x)=resinv(nep,λ=λ0,v=ones(n),logger=1,tol=1e-16);
Precomputing linsolver
iter 1 err:0.003863455199119409 λ=-1.2 + 0.0im
iter 2 err:0.0012874946992780317 λ=-1.175478914232863 + 0.0im
iter 3 err:0.016919177734890205 λ=-0.9045032212171923 + 0.0im
iter 4 err:7.884718326366283e-5 λ=-1.1998094425367551 + 0.0im
iter 5 err:9.586775788595438e-5 λ=-1.1974656536889654 + 0.0im
iter 6 err:3.7870317997772634e-5 λ=-1.1994205846695276 + 0.0im
iter 7 err:3.0752556063829646e-5 λ=-1.1985663733901704 + 0.0im
...
iter 69 err:5.647660764089489e-16 λ=-1.1989892137958458 + 0.0im
iter 70 err:5.055053391642266e-16 λ=-1.1989892137958578 + 0.0im
iter 71 err:1.781530900695102e-17 λ=-1.1989892137958473 + 0.0im
```
We will carry out some timing experiments, so let's
use the `BenchmarkTools`-package and swith off printouts in
the NEP-solver:
```julia-repl
julia> using BenchmarkTools
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),tol=1e-16);
  6.151 ms (35888 allocations: 10.65 MiB)
```
The linear system that has to be solved in every iteration
in `resinv` has a constant system matrix, and
therefore a prefactorization (typically an LU-factorization)
is useful. This is done with the `FactorizeLinSolveCreator`,
which is actually the default behaviour, so we
get no substantial difference when we specify
a creator if the type `FactorizeLinSolverCreator`.
```julia
julia> creator=FactorizeLinSolverCreator();
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator);
  6.026 ms (35879 allocations: 10.65 MiB)
```
If we do not want to use a prefactorization, you can specify
`BackslashLinSolverCreator` as your creator object.
```julia
julia> creator=BackslashLinSolverCreator();
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator);
  12.386 ms (39470 allocations: 20.29 MiB)
```
This does not use a prefactorization using an LU-factorization
and is therefore slower.


The above approach corresponded to direct methods for linear systems.
You can also use iterative methods, e.g., the GMRES-method.
The GMRES-method is available in the `GMRESLinSolverCreator`
function. Iterative methods in general need preconditioners.
We continue the example and use a diagonal preconditioner:
```julia
julia> D0=(Diagonal(compute_Mder(nep,λ0))); # Preconditioner
julia> creator=GMRESLinSolverCreator(Pl=D0, tol=1e-10);
```
All the keyword arguments in the call `GMRESLinSolverCreator`
are passed to
[`gmres!`](https://juliamath.github.io/IterativeSolvers.jl/dev/linear_systems/gmres/).
Hence, the `tol` here  specifies a termination criteria for the GMRES-method,
and `Pl` specifies the left preconditioner, in this case just a diagonal matrix.
```julia
julia> (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator,logger=1,tol=1e-16)
Precomputing linsolver
iter 1 err:0.003863455199119409 λ=-1.2 + 0.0im
iter 2 err:0.0012874946993129863 λ=-1.175478914232863 + 0.0im
iter 3 err:0.01691917773640539 λ=-0.9045032211965017 + 0.0im
iter 4 err:7.884718724884878e-5 λ=-1.1998094425350523 + 0.0im
iter 5 err:9.586776747525804e-5 λ=-1.1974656535293027 + 0.0im
...
iter 69 err:7.114086177493848e-16 λ=-1.1989892137958427 + 0.0im
iter 70 err:3.434521011774537e-16 λ=-1.1989892137958578 + 0.0im
iter 71 err:1.1060298018746236e-16 λ=-1.1989892137958504 + 0.0im
iter 72 err:1.1903842077731908e-16 λ=-1.1989892137958529 + 0.0im
iter 73 err:3.82409101411086e-16 λ=-1.1989892137958553 + 0.0im
iter 74 err:2.8385047264009487e-16 λ=-1.1989892137958476 + 0.0im
iter 75 err:2.578560246193603e-16 λ=-1.1989892137958535 + 0.0im
iter 76 err:3.8680157882039993e-16 λ=-1.198989213795848 + 0.0im
iter 77 err:1.4633231742581634e-16 λ=-1.1989892137958562 + 0.0im
iter 78 err:2.324747046412835e-16 λ=-1.198989213795853 + 0.0im
iter 79 err:1.7053137640282208e-16 λ=-1.1989892137958487 + 0.0im
iter 80 err:5.364097787291277e-17 λ=-1.1989892137958522 + 0.0im
```
The printout reveals that we need more iterations, than with a
direct method. In terms of computation
time, this approach can however still be competitive
and faster than a direct approach:
```julia-repl
julia> creator=GMRESLinSolverCreator(Pl=D0, tol=1e-2);
julia> @btime (λ,x)=resinv(nep,λ=λ0,v=ones(n),maxit=100,linsolvercreator=creator,tol=1e-16)
  8.680 ms (63237 allocations: 21.30 MiB)
```

## Your own linear solver

There are many ways to solve linear systems in Julia, e.g.,
by using package such as [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl),
[Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl)
or [KrylovMethods.jl](https://github.com/JuliaInv/KrylovMethods.jl).
These are not natively supported by NEP-PACK,
but due to the extendability of the `LinSolverCreator`-objects
specified above, you can still use them. We illustrate
the use of the implementation of BiCG in KrylovKit.


Julia's efficiency stems partially from the extensive use of types.
We need to defined new types to define our own
linear solver.

```julia
julia> struct MyLinSolverCreator <: LinSolverCreator; end
julia> struct MyLinSolver <: LinSolver;
  λ
end
```
NEP-solvers call the function `create_linsolver(linsolver,nep,λ)`,
which should return a linear solver.
This is to allow
precomputation. In our case we don't need precomputation
so we just return an instance of `MyLinSolver`.
```julia
julia> import NonlinearEigenproblems.create_linsolver # Needed since we want overload it
julia> function create_linsolver(::MyLinSolverCreator,nep,λ)
   return MyLinSolver(λ);
end
create_linsolver (generic function with 4 methods)

julia> nep=nep_gallery("dep0",50);
```
The implementation of the solver goes in the function `lin_solve`:
```julia
julia> function schur_complement_lin_solve(AA,b,n0)
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
julia> import NonlinearEigenproblems.LinSolvers.lin_solve
julia> function lin_solve(solver::MyLinSolver,b::Vector;tol=eps())
   n0=10;
   return schur_complement_lin_solve(compute_Mder(nep,solver.λ),b,n0)
end
```
You can now solve the problem by handing `MyLinSolverCreator()` to a
NEP-solver, e.g., `augnewton`:
```julia
julia> augnewton(nep,v=ones(size(nep,1)),logger=1);
iter 1 err:0.09618148118463332 λ=0.0 + 0.0im
iter 2 err:0.041321766723793736 λ=0.8989908878138995 + 0.0im
iter 3 err:0.026319557940841405 λ=1.444824857193296 + 0.0im
iter 4 err:0.004493445102203272 λ=0.9306146373776659 + 0.0im
iter 5 err:9.03013896973475e-5 λ=0.9176939995973298 + 0.0im
iter 6 err:9.801547301144594e-7 λ=0.9205599293430559 + 0.0im
iter 7 err:1.212523395823761e-10 λ=0.9205883567517947 + 0.0im
iter 8 err:6.279907189232758e-17 λ=0.9205883602865768 + 0.0im
```