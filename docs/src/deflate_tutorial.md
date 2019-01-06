# Tutorial: Computing several solutions with deflation

## Background
Several algorithms for NEPs compute one solution to the NEP
given a starting value. In many applications several
solutions are of interest. There is a trivial partial "work-around": You can try to
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
Note that starting with `λ=0` and `λ=-1` lead to the same solution.
Other solution methods do not suffer from this, e.g.,
[block Newton method](methods.md#NonlinearEigenproblems.NEPSolver.blocknewton),
[the infinite Arnoldi method](methods.md#NonlinearEigenproblems.NEPSolver.iar)
and
[nleigs](methods.md#NonlinearEigenproblems.NEPSolver.nleigs)
since they compute several solutions at once.
Another remedy is the use of *deflation*.

## Deflation in NEP-PACK

The term deflation is referring to making
something smaller (in the sense of opposite of inflating a balloon). In this case we can make the solution set smaller. We compute a solution and subsequently compute
reconstruct a deflated problem, which has the same solution as the original
problem except of the solution we have already computed.

A general deflation technique is available in NEP-PACK based on increasing
the problem size. The technique is inspired by what is described in the [PhD thesis
of Cedric Effenberger](http://sma.epfl.ch/~anchpcommon/students/effenberger.pdf).
It is implemented in the method [effenberger_deflation](transformations.md#NonlinearEigenproblems.NEPTypes.effenberger_deflation).

```julia
julia> # first compute a solution
julia> (λ1,v1)=quasinewton(nep,λ=0,v=ones(size(nep,1)))
julia> # Construct a deflated NEP where we remove (λ1,v1)
julia> dnep=effenberger_deflation(nep,λ1,v1)
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

The above procedure can be repeated by calling `effenberger_deflation` on
the deflated NEP. The procedure can be repeated such that many eigenvalues
can be computed. We now illustrate a slightly more robust variant of that
approach by constructing a deflated NEP from the original NEP.
This requires slightly more manipulation/understanding of invariant pairs which
need to be extracted for every computed solution.


```julia
function multiple_deflation(nep,λ0,p)
   n=size(nep,1);
   dnep=nep;
   S0=zeros(ComplexF64,p,p);
   V0=zeros(ComplexF64,size(nep,1),p);
   S=view(S0,1:0,1:0);
   V=view(V0,1:n,1:0);
   for k=1:p
      # Compute one solution of the deflated problem
      (λ2,v2)=quasinewton(dnep,λ=λ0,v=ones(size(dnep,1)),maxit=1000);
      # expand the invariant pair
      S0[1:k-1,k]=v2[(n+1):end];
      S0[k,k]=λ2;
      S=view(S0,1:k,1:k)
      V0[1:n,k]=v2[1:n];
      V=view(V0,1:n,1:k);
      @show S
      # Construct the deflated problem
      dnep=effenberger_deflation(nep,S,V)
   end
   return (S,V)
end
```

We can now compute several solutions by calling `multiple_deflation`.
Note that we use the same starting eigenvalue for all eigenvalues: `0.5im`. It has
to be complex in this case, since if it was real, we would not find complex solution and this problem only has two real eigenvalues.
```julia
julia> nep=nep_gallery("dep0");
julia> (S,V)=multiple_deflation(nep,0.5im,3)
S = Complex{Float64}[-0.358719+1.33901e-14im]
S = Complex{Float64}[-0.358719+1.33901e-14im -0.769266-0.728303im; 0.0+0.0im 0.834735+1.25838e-14im]
S = Complex{Float64}[-0.358719+1.33901e-14im -0.769266-0.728303im -0.735867-0.43166im; 0.0+0.0im 0.834735+1.25838e-14im 0.570725-0.153773im; 0.0+0.0im 0.0+0.0im -0.0409352+1.48601im]
```

The matrix pair `(S,V)` is a partial Schur factorization of the NEP. This can be
seen from the fact [compute_MM](types.md#NonlinearEigenproblems.NEPCore.compute_MM) vanishes
for for `(S,V)`:
```julia
julia> norm(compute_MM(nep,S,V))
4.341002168663845e-13
```
The eigenvalues of `S` are eigenvalues of the original NEP, and we can find the eigenpairs by
diagonalizing `S`:
```julia
julia> (Λ,P)=eigen(S);
julia> VV=V*P;  # Construct the eigenvector matrix
julia> Λ # The computed eigenvalues
3-element Array{Complex{Float64},1}:
  -0.3587189459686267 + 1.339010598765711e-14im
   0.8347353572199371 + 1.2583846244197297e-14im
 -0.04093521177096655 + 1.4860115309416284im
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
