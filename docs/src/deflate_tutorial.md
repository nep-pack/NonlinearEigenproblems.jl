# Tutorial: Computing several solutions with deflation

## Background
Several algorithms for NEPs compute one solution to the NEP
given a starting value. In many applications several
solutions are of interest. Let us first consider the trivial partial
"work-around": You can try to
run an algorithm which computes one eigenvalue twice with
different starting values, e.g., quasinewton as in this
example:
```julia-repl
julia> using NonlinearEigenproblems
julia> nep=nep_gallery("dep0",4);
julia> (λ1,v1)=quasinewton(nep,λ=0,v=ones(size(nep,1)),maxit=500)
(-0.2094378352960786 + 0.0im, Complex{Float64}[2.2479093650910866 + 0.0im, 3.208352895087154 + 0.0im, 0.6628450056308428 + 0.0im, 17.946407917249605 + 0.0im])
julia> (λ2,v2)=quasinewton(nep,λ=1,v=ones(size(nep,1)),maxit=500)
(0.2966714721676867 + 0.0im, Complex{Float64}[-0.8369951877176647 + 0.0im, -5.749143077718012 + 0.0im, 5.720770822961643 + 0.0im, 10.495353352793199 + 0.0im])
```
This simple approach often suffers from the problem called *reconvergence* (we obtain the
same solution again) or solutions of interest may be missed. In this case we get
reconvergence when we use starting value `-1`:
```julia-repl
julia> (λ3,v3)=quasinewton(nep,λ=-1,v=ones(size(nep,1)),maxit=500)
(-0.20943783529618362 + 0.0im, Complex{Float64}[-0.6110572600894503 + 0.0im, -0.8721380674494782 + 0.0im, -0.18018353377334895 + 0.0im, -4.878436391013284 + 0.0im])
```
Note that applying the algorithm with starting values `λ=0` and `λ=-1` lead to the same solution.
Other solution methods do not suffer from this, e.g.,
[block Newton method](methods.md#NonlinearEigenproblems.NEPSolver.blocknewton),
[the infinite Arnoldi method](methods.md#NonlinearEigenproblems.NEPSolver.iar)
and
[nleigs](methods.md#NonlinearEigenproblems.NEPSolver.nleigs)
since they compute several solutions at once.
Another attempt to remedy reconvergence
is to use the technique called *deflation*. See also
[the manual page on deflation](deflation.md).

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


In NEP-PACK, this type of deflation is implemented in the function [`deflate_eigpair`](@ref),
which takes a NEP and an eigenpair as input and returns a new NEP.
```julia-repl
julia> # Construct a deflated NEP where we remove (λ1,v1)
julia> dnep=deflate_eigpair(nep,λ1,v1);
julia> # The dnep is a new NEP but with dimension increased by one
julia> size(nep)
(4, 4)
julia> size(dnep)
(5, 5)
```
We now illustrate that we can avoid reconvergence:
```julia
julia> (λ4,v4)=quasinewton(dnep,λ=-1,v=ones(size(dnep,1)),maxit=500)
(0.29667147216767376 + 0.0im, Complex{Float64}[-11.767671406737819 + 0.0im, -43.86197116968253 + 0.0im, 31.9938464980679 + 0.0im, 8.133682253178579 + 0.0im, -28.114795306465478 + 0.0im])
```
Note: In contrast to the initial example, starting value `λ=-1` does *not* lead to converge to the eigenvalue we obtained from starting value `λ=0`.

The computed solution is indeed a solution to the original NEP since ``M(λ4)`` is singular:
```julia
julia> minimum(svdvals(compute_Mder(nep,λ4)))
4.166120681513672e-14
```
In fact, you can even start with the first starting value `λ=0`, and get a new solution
```julia
julia> quasinewton(dnep,λ=0,v=ones(size(dnep,1)),maxit=500)
(-1.3640414062700734 + 0.0im, Complex{Float64}[-1.0976664883572566e307 + 0.0im, -2.8870394809137054e307 + 0.0im, 2.1189933442957902e307 + 0.0im, 5.753536946292879e306 + 0.0im, -1.9207807191677339e307 + 0.0im])
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
Note that we use the same starting guess, `1im`, for all eigenvalues.
```julia-repl
julia> (Λ,VV)=multiple_deflation(nep,1im,3)
(Complex{Float64}[-0.20943783529608836 - 1.1876898672667347e-13im, 0.09916136114196937 + 1.5040449444406283im, 0.2966714721675515 - 1.1055942412512364e-13im], Complex{Float64}[0.08606692232135099 - 0.0868831949158175im -0.33250443251734496 + 0.6897411986346054im -0.04563655283979674 + 0.04339934493508566im; 0.12283994350000517 - 0.12400497736757769im 0.2101087999777705 - 0.18967830722177997im -0.31346783792781846 + 0.29810092957798906im; 0.025378705430330793 - 0.02561940117237799im 0.41435154273390934 + 0.36869913908327107im 0.31192086140711456 - 0.2966297893745159im; 0.6871238316577664 - 0.6936406250750684im -0.15130500943732544 - 0.05527070830243071im 0.5722514954516473 - 0.5441984219948924im])
```
The values in `Λ` and `VV` are eigenpairs:
```julia
julia> norm(compute_Mlincomb(nep,Λ[1],VV[:,1]))
1.7819713566771836e-13
julia> norm(compute_Mlincomb(nep,Λ[2],VV[:,2]))
1.2888961114892419e-13
julia> norm(compute_Mlincomb(nep,Λ[3],VV[:,3]))
1.5058274131661697e-13
```


![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_DEFLATION)
