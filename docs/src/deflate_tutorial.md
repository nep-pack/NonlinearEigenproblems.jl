# Tutorial: Computing several solutions with deflation

## Background
Several algorithms for NEPs compute one solution to the NEP
given a starting value. In many applications several
solutions are of interest. Let us first consider the trivial partial
"work-around": You can try to
run an algorithm which computes one eigenvalue twice with
different starting values, e.g., quasinewton as in this
example:
```julia
julia> nep=nep_gallery("dep0");
julia> (λ1,_)=quasinewton(nep,λ=0,v=ones(size(nep,1)))
(-0.3587189459686377 + 0.0im, Complex{Float64}[4.41411+0.0im, -2.22171+0.0im, 4.31544+0.0im, -7.76501+0.0im, -9.51261+0.0im])
julia> (λ2,_)=quasinewton(nep,λ=1im,v=ones(size(nep,1)))
(-0.04093521177097875 + 1.4860115309416284im, Complex{Float64}[-3.28271+11.7399im, 5.08623-8.05479im, 7.16697-6.25547im, -2.69349+4.63954im, -9.91065+14.4678im])
```
This simple approach often suffers from the problem called *reconvergence* (we obtain the
same solution again) or solutions of interest may be missed. In this case we get
reconvergence when we use starting value `-1`:
```julia
julia> (λ3,_)=quasinewton(nep,λ=-1,v=ones(size(nep,1)))
(-0.358718945968621 + 0.0im, Complex{Float64}[-6.65881+0.0im, 3.35151+0.0im, -6.50997+0.0im, 11.7137+0.0im, 14.3501+0.0im])
```
Note that applying the algorithm with starting values `λ=0` and `λ=-1` lead to the same solution.
Other solution methods do not suffer from this, e.g.,
[block Newton method](methods.md#NonlinearEigenproblems.NEPSolver.blocknewton),
[the infinite Arnoldi method](methods.md#NonlinearEigenproblems.NEPSolver.iar)
and
[nleigs](methods.md#NonlinearEigenproblems.NEPSolver.nleigs)
since they compute several solutions at once.
Another remedy is the use of *deflation*.

## Deflation in NEP-PACK

The term deflation is referring to making
something smaller (in the sense of opposite of inflating a balloon). In this case we can make the solution set smaller. We compute a solution and subsequently
construct a deflated problem, which has the same solutions as the original
problem except of the solution we have already computed.

A general solver independent deflation technique is available in NEP-PACK based on increasing
the problem size.
(There are also NEP-solver deflation techniques incoprorated in, e.g., in [the nonlinear Arnoldi method](methods.md#NonlinearEigenproblems.NEPSolver.nlar) and [the Jacobi-Davidson method](methods.md#NonlinearEigenproblems.NEPSolver.jd_betcke).)
The solver independent technique is inspired by what is described in the [PhD thesis
of Cedric Effenberger](http://sma.epfl.ch/~anchpcommon/students/effenberger.pdf).
It is implemented in the method [effenberger_deflation](transformations.md#NonlinearEigenproblems.NEPTypes.effenberger_deflation).

```julia
julia> # first compute a solution
julia> (λ1,v1)=quasinewton(nep,λ=0,v=ones(size(nep,1)))
julia> # Construct a deflated NEP where we remove (λ1,v1)
julia> dnep=deflate_eigpair(nep,λ1,v1)
julia> # The dnep is a new NEP but with dimension increased by one
julia> size(nep)
(5, 5)
julia> size(dnep)
(6, 6)
```
We now illustrate that we can avoid reconvergence:
```julia
julia> (λ4,v4)=quasinewton(dnep,λ=-1,v=ones(size(dnep,1)),maxit=1000)
(0.8347353572199264 + 0.0im, Complex{Float64}[10.6614+0.0im, 0.351814+0.0im, -0.940539+0.0im, 1.10798+0.0im, 3.53392+0.0im, -0.447213+0.0im])
```
Note: In contrast to the initial example, starting value `λ=-1` does *not* lead to converge to the eigenvalue we obtained from starting value `λ=0`.

The computed solution is indeed a solution to the original NEP since `M(λ4)` is singular:
```julia
julia> minimum(svdvals(compute_Mder(nep,λ4)))
1.2941045763733582e-14
```
In fact, you can even start with the first starting value `λ=0`, and get a new solution
```julia
julia> quasinewton(dnep,λ=0,v=ones(size(dnep,1)),maxit=1000)
(0.8347353572199577 + 0.0im, Complex{Float64}[9.28596+0.0im, 0.306425+0.0im, -0.819196+0.0im, 0.965031+0.0im, 3.07799+0.0im, -0.389516+0.0im])
```

## Repeated deflation

The above procedure can be repeated by calling `deflate_eigpair` on
the deflated NEP. This effectively deflates another eigenpair
(but without creating a recursive deflated nep structure).


```julia
function multiple_deflation(nep,λ0,p)
   n=size(nep,1);
   dnep=nep;
   for k=1:p
      # Compute one solution of the deflated problem
      (λ2,v2)=quasinewton(dnep,λ=λ0,v=ones(size(dnep,1)),maxit=1000);
      # expand the invariant pair
      dnep=deflate_eigpair(dnep,λ2,v2)
   end
   return get_deflated_eigpairs(dnep);

end
```

We can now compute several solutions by calling `multiple_deflation`.
Note that we use the same starting eigenvalue for all eigenvalues: `0.5im`. It has
to be complex in this case, since if it was real, we would not find complex solution and this problem only has two real eigenvalues.
```julia
julia> nep=nep_gallery("dep0");
julia> (Λ,VV)=multiple_deflation(nep,0.5im,3)
(Complex{Float64}[-0.358719+1.33901e-14im, 0.834735+7.05729e-15im, -0.0409352+1.48601im], Complex{Float64}[-0.0148325-0.316707im -0.670282+0.268543im -0.41261+0.229832im; 0.00746549+0.159405im -0.0881321+0.0353094im 0.360381-0.0796982im; … ; 0.0260924+0.557131im -0.298976+0.119782im -0.201138+0.0524051im; 0.0319648+0.68252im -0.528234+0.211633im -0.668441+0.121828im])
```
The values in `Λ` and `VV` are eigenpairs:
```julia
julia> norm(compute_Mlincomb(nep,Λ[1],VV[:,1]))
2.0521012310648373e-13
julia> norm(compute_Mlincomb(nep,Λ[2],VV[:,2]))
2.8707903010898464e-13
julia> norm(compute_Mlincomb(nep,Λ[3],VV[:,3]))
1.883394132275381e-13
```

## The theory in the background

The deflation is based on a theory for NEP essentially stating that
if ``(s,x)`` is an eigenpair, then the extended nonlinear eigenvalue problem
```math
T(λ):=\begin{bmatrix}M(λ)&M(λ)x(s-λ)^{-1}\\ x^T & 0\end{bmatrix}
```
has the same eigenvalues as the original problem (under certain conditions
quite general conditions which are assumed to be satisfied). More
eigenpairs can be deflated with techniques of partial Schur
factorizations (which the user does not need to use). When we create
a deflated NEP, we create the NEP `T`.

There are several ways to represent the ``T``, which is why deflation has several
modes. If you run
```julia
julia> dnep=deflate_eigpair(nep,λ1,v1,mode=:SPMF)
```
the `dnep` will be of the type `AbstractSPMF`. More precisely, if
```math
M(λ)=A_1f_1(λ)+\cdots+A_mf_m(λ)
```
the deflated NEP will be
```math
T(λ)=
\begin{bmatrix}A_1&0\\0 & 0\end{bmatrix}f_1(λ)+\cdots+
\begin{bmatrix}A_m&0\\0 & 0\end{bmatrix}f_m(λ)+
\begin{bmatrix}0&A_1x\\0 & 0\end{bmatrix}\frac{f_1(λ)}{s-λ}+\cdots+
\begin{bmatrix}0&A_mx\\0 & 0\end{bmatrix}\frac{f_m(λ)}{s-λ}+
\begin{bmatrix}0&0\\x^T & 0\end{bmatrix}
```
Clearly, the deflated NEP will have more SPMF-terms, and
the `mode=:SPMF`, is not recommended if you have many SPMF-terms.
(Some additional exploitation is however implemented, since we can use
the fact that the introduced terms are of low rank, and
therefore naturally represented as a `LowRankFactorizedNEP`.)

If you select `mode=:Generic`, the compute functions are implemented
without the use of SPMF, and can be more efficient
if the NEP has many SPMF-terms.
When `mode=:MM` the compute-functions are all implemented
by calls to `compute_MM`. This will not be efficient if
`compute_Mder(nep,λ,der)` where  `der>0` is needed.

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_DEFLATION)
