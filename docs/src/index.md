

# NEP-PACK

NEP-PACK is a package with implementations of methods to solve and to manipulate
nonlinear eigenvalue problems of
the type: Find ``(λ,v)\in\mathbb{C}\times\mathbb{C}^n`` such that
```math
M(λ)v=0
```
and ``v\neq 0``.


## Getting started

Install it as a registered  package in Julia's REPL package mode by
typing `] add Nonline...`:
```
julia> ]
(v1.0) pkg> add NonlinearEigenproblems
```
Then we can start to load the NEP-PACK package
```julia-repl
julia> using NonlinearEigenproblems
```
As a first example we will solve the NEP associated with the matrix polynomial

```math
M(λ)=\begin{bmatrix}1&3\newline5&6\end{bmatrix}+
λ\begin{bmatrix}3&4\newline6&6\end{bmatrix}+
λ^2\begin{bmatrix}1&0\newline0&1\end{bmatrix}
```
The following code creates this NEP, by constructing an object called
[`PEP`](@ref), an abbreviation for polynomial eigenvalue problem.
It subsequently solves it using the NEP solution method implemented
in the NEP-solver [`polyeig`](@ref):
```julia-repl
julia> A0=[1.0 3; 5 6]; A1=[3.0 4; 6 6]; A2=[1.0 0; 0 1.0];
julia> nep=PEP([A0,A1,A2])
PEP(2, Array{Float64,2}[[1.0 3.0; 5.0 6.0], [3.0 4.0; 6.0 6.0], [1.0 0.0; 0.0 1.0]])
julia> λ,v=polyeig(nep)
(Complex{Float64}[1.36267+0.0im, -0.824084+0.280682im, -0.824084-0.280682im, -8.7145+0.0im], Complex{Float64}[-1.0+0.0im 0.739183-0.196401im 0.739183+0.196401im 0.627138+0.0im; 0.821812+0.0im -0.501408-0.375337im -0.501408+0.375337im 1.0+0.0im])
```
You have now solved your first nonlinear eigenvalue problem with NEP-PACK.

In order to verify that we have a solution, we can check that  ``M(λ)`` is singular,
with a singular vector ``v`` such that ``M(λ)v=0``:
```julia-repl
julia> λ1=λ[1]; v1=v[:,1];
julia> using LinearAlgebra # the norm-function is in this Julia package
julia> norm(A0*v1+λ1*A1*v1+λ1^2*v1)/norm(v1)
1.1502634749464687e-14
```

!!! tip
    MATLAB users: Do you have a NEP defined in MATLAB? You can solve MATLAB-defined NEPs with this package.  See [the MATLAB tutorial](tutorial_matlab1.md). We also have some MATLAB implementations of the solvers in NEP-PACK in a [separate repository](https://github.com/nep-pack/NEP-PACK-matlab-reference).


## Accessing more complicated applications

We have made benchmark examples available through the function [`nep_gallery`](gallery.md#NonlinearEigenproblems.nep_gallery):

```julia-repl
julia> nep=nep_gallery("dep0",100);
julia> size(nep)
(100, 100)
julia> λ,v=mslp(nep,tol=1e-10);
julia> λ
0.05046248970129549 - 7.60684247532422e-16im
julia> size(v)
(100,)
julia> resnorm=norm(compute_Mlincomb(nep,λ,v))
5.178780131881974e-13
```
Information about the gallery can be found by typing `?nep_gallery`.
The second argument in the call to `nep_gallery` is a problem parameter,
in this case specifying that the  size of the problem should be `100`.
The example solves the problem with the NEP-algorithm [`MSLP`](methods.md#NonlinearEigenproblems.NEPSolver.mslp).
The parameter `tol` specifies the
tolerance for iteration termination.

!!! note
    All the NEP-solvers have considerble documentation easily available.
    Every NEP-solver has documentation accompanied with at least one example,
    and references to corresponding research papers, which we strongly recommend you
    to cite if you use the method.
    This is available to you in Julia's repl-prompt. Type `?mslp` and you will see
    an example how to use `mslp` and that citation credit should go to *A. Ruhe,
    Algorithms for the nonlinear eigenvalue problem, SIAM J. Numer. Anal.
    10 (1973) 674-689*. This documentation is the same as the online documentation
    under the tab [NEP-solvers](methods.md).



## A model of a neuron

The following (delay) differential equation models the interaction of two neurons
```math
\dot{x}_1(t)=-\kappa x_1(t)+\beta\tanh(x_1(t-\tau_3))+a_1\tanh(x_2(t-\tau_2))
```
```math
\dot{x}_2(t)=-\kappa x_2(t)+\beta\tanh(x_2(t-\tau_3))+a_2\tanh(x_1(t-\tau_1))
```
See [L. P. Shayer and S. A. Campbell.  Stability, bifurcation and multistability in a system of two
coupled neurons with multiple time delays. SIAM J. Applied Mathematics , 61(2):673–700, 2000](https://www.jstor.org/stable/3061744?seq=1#page_scan_tab_contents). It is
also available as a first demo in [DDE-BIFTOOL](https://sourceforge.net/projects/ddebiftool/).
The linear stability analysis of this problem requires the solution
of a nonlinear eigenvalue problem
```math
M(λ)=-λI+A_0+A_1e^{-\tau_1λ}+A_2e^{-\tau_2λ}+A_3e^{-\tau_3λ}
```
where the matrices are the Jacobian at the stationary solution.
For the zero stationary solution, the matrices are
```julia-repl
kappa=0.5; a2=2.34; a1=1; beta=-1;
A0=-kappa*[1 0; 0 1];
A1=a2*[0 0; 1 0];
A2=a1*[0 1; 0 0];
A3=beta*[1 0; 0 1];
```
We can now create the nonlinear eigenvalue problem and determine the stability
by first creating the problem
```julia-repl
julia> tauv=[0;0.2;0.2;1.5];
julia> dep=DEP([A0, A1,   A2, A3],tauv);
```
The constructor  [`DEP`](@ref) is an abbreviation for a delay eigenvalue problem, which
is a NEP with exponential terms stemming from the stability
analysis of a delay-differential equation. See [Types and data-structures](types.md) for other NEP-types.
You can now solve this NEP, for instance,
with the [infinite Arnoldi method](methods.md#NonlinearEigenproblems.NEPSolver.iar_chebyshev):
```julia-repl
julia> λ,V=iar_chebyshev(dep,maxit=100); # This takes some time the first time is run due to JIT-compiler
```
The figure in a demo of DDE-BIFTOOL <http://ddebiftool.sourceforge.net/demos/neuron/html/demo1_stst.html#3> can be directly generated by


```@setup neuron
using Gadfly
set_default_plot_size(12cm, 12cm)
```

```@example neuron
using Gadfly
λ=[ -0.09712795241565722 + 2.612885243197631e-19im # hide
         0.30886599775839135 + 4.146563548756125e-18im # hide
        -0.45584765486526174 + 1.6884551234089458im # hide
         -0.4558476548652613 - 1.6884551234089418im # hide
         -0.8832708076887316 + 5.325050575287575im # hide
         -0.8832708076887288 - 5.3250505752875625im] # hide
plot(x=real.(λ),y=imag.(λ), Guide.xlabel("real(λ)"), Guide.ylabel("imag(λ)"))
```

!!! tip
    This problem is also available in the `Gallery` by calling `dep=nep_gallery("neuron0")`. Most of the NEPs constructed in the tutorials are also available in corresponding gallery problems. See all gallery problems under [NEP Gallery](gallery.md). In particular, note that the problems in the Berlin-Manchester collection of problems NLEVP are also [directly available](gallery.md#Berlin-Manchester-collection-1).

## The "gun" benchmark problem

One of the most common benchmark problems for NEPs is the so-called "gun"-problem.
It models an electromagnetic cavity, and it is directly available in the NEP-PACK
gallery.
(See [gallery](gallery.md#NonlinearEigenproblems.Gallery.nep_gallery) references or type `?nep_gallery` at the repl-prompt.) This is how you can set it up and solve it with the [block Newton method](methods.md#NonlinearEigenproblems.NEPSolver.blocknewton):

```julia-repl
julia> nep=nep_gallery("nlevp_native_gun");
julia> n=size(nep,1)
julia> S=150^2*[1.0 0; 0 1]; V=[[1 0; 0 1]; zeros(n-2,2)];
julia> (Z,X)=blocknewton(nep,S=S,X=V,logger=1,armijo_factor=0.5,maxit=20)
Iteration 1: Error: 6.081316e+03
Iteration 2: Error: 1.701970e-02 Armijo scaling=0.031250
Iteration 3: Error: 1.814887e-02 Armijo scaling=0.250000
...
Iteration 13: Error: 6.257442e-09
Iteration 14: Error: 2.525942e-15
```
This algorithm returns a partial Schur factorization
of the NEP, and therefore the eigenvalues of the small matrix
`Z` are eigenvalues of the problem. An eigenpair of the NEP
can be extracted by diagonalizing:
```julia-repl
julia> using LinearAlgebra
julia> (Λ,P)=eigen(Z);
julia> VV=X*P;  # Construct the eigenvector matrix
julia> v=VV[:,1]; λ=Λ[1]
61330.208714730004 + 63185.15983933589im
julia> norm(compute_Mlincomb(nep,λ,v)) # Very small residual
1.8270553408452648e-16
```

If you use the NEP-algorithms for research, please
give the author of the algorithm credit by citiation. The
recommended citation can be found in the function
documentation, e.g., `?blocknewton`.

## Your own NEP nonlinearity

As an application researcher, we recommend that you first try to
express your problem in the following form since it
gives access to several efficient routines associated with the NEP,
in turn making it possible to use many NEP-solvers. A problem that can be expressed as a (short) **S** um of **P** roducts of **M** atrices and **F** unctions
can be represented with the objects of type [`SPMF_NEP`](@ref)
in NEP-PACK. For instance, a problem with three terms
```math
M(λ) = A+λB+e^{\sin(λ/2)}C
```
can be created by
```julia-repl
julia> A=(1:4)*(1:4)'+I; B=diagm(1 => [1,2,3]); C=ones(4,4);
julia> f1= λ-> one(λ);
julia> f2= λ-> λ;
julia> f3= λ-> exp(sin(λ/2));
julia> nep=SPMF_NEP([A,B,C],[f1,f2,f3]);
```
The NEP is solved by using the NEP-object as a parameter in a call to an algorithm, e.g.,
```julia-repl
julia> v0 = 0.1*[1,-1,1,-1];
julia> λ,v=quasinewton(nep,λ=4,v=v0)
(3.1760990071435193 + 0.0im, Complex{Float64}[2.892363187499394 + 0.0im, -1.6573097795628646 + 0.0im, 0.00729776922332883 + 0.0im, -0.09002519738673213 + 0.0im])
```
As usual, you can check that we computed a sensible solution:
```julia-repl
julia> (A+B*λ+C*exp(sin(λ/2)))*v
4-element Array{Complex{Float64},1}:
  -3.489601657766542e-12 + 0.0im
 -1.0118303586944344e-12 + 0.0im
  -9.480334553029193e-13 + 0.0im
  -5.912084880273861e-13 + 0.0im
```
!!! note
    The functions `f1`,`f2` and `f3` in the example above have to be defined for scalar values and for matrices (in the [matrix function](https://en.wikipedia.org/wiki/Matrix_function) sense, not elementwise sense). This is the reason `f1` needs to be defined as `one(λ)`, instead of just `1`. Fortunately, many elementary functions in Julia already have matrix function implementations, e.g., `exp([1 2 ; 3 4])` will return the matrix exponential of the given matrix.



## Chebyshev interpolation

In applications, NEP-nonlinearities may be complicated to implement.
Directly using the SPMF-functionality where every function needs to be defined in a
matrix function sense may require too much work.
In this case you may want to use an approximation method to
create a new different NEP object for which the matrix functions are
easy to implement (or directly available in the package).
We illustrate this property with NEP-PACKs Chebyshev interpolation feature.


Suppose you have the following NEP, which requires a Bessel function.
The Bessel function is analytic, but its matrix function is not
easily available.
```julia-repl
julia> using SpecialFunctions; # for the besselj
julia> fv=Vector{Function}(undef,m);
julia> Av=Vector{Matrix{Float64}}(undef,3)
julia> fv[1]=s->one(s);
julia> Av[1]=[ -2.0  -1.0   8.0; -1.0  0  -1.0;   -2.0   -1.0  -2.0];
julia> fv[2]=s->s;
julia> Av[2]=[4.0 -7.0  14.0; 8.0  9.0 -13.0; -1.0 -1.0    10.0];
julia> fv[3]=s->besselj(0, s);
julia> Av[3]=[-7.0 -0.0 -9.0; 8.0  3.0 -3.0;  0.0 13.0  2.0]
```
We use `SPMF_NEP` again, but in order to suppress a warning message
indicating that evaluation with a  matrix function
is not available we use the keyword `check_consistency=false`.
```julia-repl
julia> nep=SPMF_NEP(Av,fv,check_consistency=false);
```

Note that we cannot directly use the `nep` object with most NEP-solvers, since
we did not provide a matrix function implementation for `besselj`. Any method
requiring a derivative will just throw an error message that a matrix function is not defined.
Let us now construct an interpolating Chebyshev polynomial,
which we can use instead (since its matrix functions are trivial).
The command [`ChebPEP`](@ref), by default  interpolates a NEP in
the interval `[-1,1]` using Chebyshev points and represent the approximation in a Chebyshev basis:
```julia-repl
julia> cheb=ChebPEP(nep,9,cosine_formula_cutoff=9);
```

We can now use an arbitrary method to try to solve this problem, e.g.,
the [`newtonqr`](@ref) method.
```julia-repl
julia> (λ,v)=newtonqr(cheb,λ=0.0,logger=1)
iter 1 err:0.20552458291903797 λ=0.0 + 0.0im
iter 2 err:0.10317368136012978 λ=-3.1180031985377803 + 0.0im
iter 3 err:0.03898166871714645 λ=-0.5814386400379581 + 0.0im
iter 4 err:0.001421286693333467 λ=-0.4572312118506711 + 0.0im
iter 5 err:1.599526685190599e-6 λ=-0.46101438033594805 + 0.0im
iter 6 err:1.9383172515233692e-12 λ=-0.4610101105535983 + 0.0im
iter 7 err:2.1034235144362163e-17 λ=-0.4610101105484241 + 0.0im
(-0.4610101105484241 + 0.0im, Complex{Float64}[-0.597958+0.0im, 0.322148+0.0im, 1.0+0.0im], Complex{Float64}[-0.257712+0.0im, -0.964465+0.0im, -0.0582387+0.0im])
```
This solved the interpolated problem quite accurately,
which turns out to be a reasonable approximation
of the original problem:
```julia-repl
julia> norm(compute_Mlincomb(nep,λ,v))
1.148749763351579e-9
```
The function `compute_Mlincomb` returns the evaluation of `M(λ)*v`; see
[the manual section for compute functions](compute_functions.md).

## What now?

Now you are ready to try out
one of our tutorials
[on artificial boundary conditions](movebc_tutorial.md),
[boundary element method](bemtutorial.md),
[contour integration](tutorial_contour.md),
or
[deflation](deflate_tutorial.md).
See also the other tutorials (in the side-bar),
or have a look at the examples
in [NEP-solvers](methods.md) and  [NEP Gallery](gallery.md).


## How do I cite it?


We have a [preprint for this work](https://arxiv.org/abs/1811.09592). If you find this software useful please cite this preprint by using this citation data:
```bibtex
@Misc{,
  author = 	 {E. Jarlebring and M. Bennedich and G. Mele and E. Ringh and P. Upadhyaya},
  title = 	 {{NEP-PACK}: A {Julia} package for nonlinear eigenproblems},
  year = 	 {2018},
  note = 	 {https://github.com/nep-pack},
  eprint = 	 {arXiv:1811.09592},
}
```
If you use a specific NEP-solver, please also give credit to the algorithm researcher.
Reference to a corresponding algorithm paper can be found by in, e.g., by writing `?resinv`.


![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC)
