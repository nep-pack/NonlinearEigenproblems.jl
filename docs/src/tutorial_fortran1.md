# Tutorial: Solving a NEP defined in fortran

## A problem defined in fortran

A situation may arise where you  have to
(or have the opportunity to) work with fortran code.
This is not as uncommon as many think, mostly
due to the legacy software in many engineering disciplines.
The Julia language is designed
with interoperability in mind. Don't let some
fortran code scare you.
The following tutorial illustrates interoperability
in Julia and how to use it in NEP-PACK.

!!! note
    To work with NEPs defined in fortran you need to compile your fortran code.
    This tutorial is written for Ubuntu Linux and
    [GNU fortran](https://gcc.gnu.org/wiki/GFortran).
    However, it is possible also under the Windows OS.

We assume our NEP is defined in fortran code and
defines the problem
```math
M(\lambda)=A_0+\lambda^3e_ne_1^T-\exp(\lambda)e_1e_n^T.
```
where ``A_0`` is a finite difference approximation of a scaled
Laplacian matrix. The problem can be naturally represented
in sparse format, which we will also take advantage of.

The fortran implementation of the problem
is given in the following subroutine which
computes three vectors `I`, `J` and `F`, where `I` and `J`
correspond to row and column pointers and `F` the value
of the sparse matrix.
The variable `s` is ``λ``, i.e., the evaluation point.
The input `der` determines which derivative
of ``M(λ)`` should be computed.

!!! tip
    Later in this tutorial we look at what can be done if the
    [derivatives are not easily available](tutorial_fortran1.md#Implementation-in-NEP-PACK:-basic-usage,-no-derivatives-1).

This is the implementation which we put in `myproblem.f95`:
```fortran
subroutine mder(s,n,der,I,J,F)
  real*8, intent(in) :: s
  integer*8, intent(in) :: n
  integer*8, intent(in) :: der
  integer*8, intent(out), dimension(3*n):: I
  integer*8, intent(out), dimension(3*n):: J
  real*8, intent(out), dimension(3*n):: F
  integer*8 :: p
  real*8 :: factor
  if (der==0) then
     factor=1;
  else
     factor=0;
  end if
  do p = 1, n
     I(p) = p
     J(p) = p
     F(p) = 2.0*factor;
  end do
  do p = 1, n-1
     I(n+p) = p
     J(n+p) = p+1
     F(n+p) = -1.0*factor;
     I(2*n+p) = p+1
     J(2*n+p) = p
     F(2*n+p) = -1.0*factor;
  end do
  I(2*n)=n;
  J(2*n)=1;
  if (der == 0) then
     F(2*n)=s*s*s;
  else if (der == 1) then
     F(2*n)=3*s*s;
  else if (der == 2) then
     F(2*n)=3*2*s;
  else if (der == 3) then
     F(2*n)=3*2;
  else
     F(2*n)=0;
  end if
  I(3*n)=1;
  J(3*n)=n;
  F(3*n)=-exp(s);
end subroutine mder
```

## Compile and call the code


Compile the code to a shared object file with the command `gfortran`:
```bash
$ gfortran -shared -fPIC -o myproblem.so myproblem.f95
```
(Under the Windows OS, you would want to compile the code to a dll-file.)
In Julia, you can now call this routine using the `Libdl`
package:
```julia
using NonlinearEigenproblems, Libdl;
mylib=Libdl.dlopen("./myproblem.so");
λ=0.3;
der=0;
n=3; # Problem size
I=Vector{Int}(undef,3*n); # 3*n nnz elements in matrix
J=Vector{Int}(undef,3*n); # 3*n nnz elements in matrix
F=Vector{Float64}(undef,3*n); # 3*n nnz elements in matrix
# This is the call to the fortran code
# Note that :mder_ is a reference to a fortran subroutine:
# it must be lower-case and  a _ should be appended
ccall(Libdl.dlsym(mylib,:mder_), Nothing,
   (Ref{Float64}, Ref{Int},Ref{Int},  Ptr{Int}, Ptr{Int}, Ptr{Float64}),
   λ, n, der, I, J, F)
```
The above code sets vectors `I`, `J` and `F` such that they
represent a sparse matrix. The sparse matrix can be
constructed with the `sparse` command:
```julia-repl
julia> using SparseArrays
julia> A=sparse(I,J,F)
3×3 SparseMatrixCSC{Float64,Int64} with 9 stored entries:
  [1, 1]  =  2.0
  [2, 1]  =  -1.0
  [3, 1]  =  0.027
  [1, 2]  =  -1.0
  [2, 2]  =  2.0
  [3, 2]  =  -1.0
  [1, 3]  =  -1.34986
  [2, 3]  =  -1.0
  [3, 3]  =  2.0
julia> Matrix(A)
3×3 Array{Float64,2}:
  2.0    -1.0  -1.34986
 -1.0     2.0  -1.0
  0.027  -1.0   2.0
```

## Implementation in NEP-PACK: basic usage
We saw above how to compute a derivative matrix with a fortran call.
This is sufficient to define a NEP-object in NEP-PACK using
the `Mder_NEP` type.

```julia
n=100;
# A function which allocates vectors and calls fortran,
# and returns a sparse matrix
function my_Mder(λ::Float64,der::Int=0)
  # Index vectors: Length 3*n since we have 3n nnz elements in matrix
  I=Vector{Int}(undef,3*n);
  J=Vector{Int}(undef,3*n);
  F=Vector{Float64}(undef,3*n);
  ccall(Libdl.dlsym(mylib,:mder_), Nothing,
     (Ref{Float64}, Ref{Int},Ref{Int},  Ptr{Int}, Ptr{Int}, Ptr{Float64}),
     λ, n, der, I, J, F)
  return sparse(I,J,F);
end
nep=Mder_NEP(n,my_Mder);
```
With the NEP defined, we can use, e.g.,
[`quasinewton`](methods.md#NonlinearEigenproblems.NEPSolver.quasinewton), to solve it:
```julia-repl
julia> quasinewton(Float64,nep,λ=-1.8,v=ones(n), logger=1);
Precomputing linsolver
Iteration:  1 errmeasure:4.903565024143569095e-01, λ=-1.8
Iteration:  2 errmeasure:8.776860766232853772e-02, λ=-1.3816406142423465
Iteration:  3 errmeasure:6.109070850428219984e-02, λ=-2.0060080798679913
...
Iteration: 11 errmeasure:5.305001776886219717e-12, λ=-1.7940561686588974
Iteration: 12 errmeasure:2.895637837297152945e-14, λ=-1.7940561686787597
Iteration: 13 errmeasure:3.874312247075750238e-16, λ=-1.7940561686786516
```


## Implementation in NEP-PACK: basic usage, no derivatives

In the above example, all the derivatives of ``M(λ)`` were
easy to compute by hand and made available in the fortran subroutine.
In many applications, the nonlinearity is not so simple,
and its derivatives may require man-hours to analyze and implement,
or may be very computationally expensive.

Most NEP-algorithms in NEP-PACK do require the derivative
(except for
certain versions of
[`nleigs`](methods.md#NonlinearEigenproblems.NEPSolver.nleigs),
[`broyden`](methods.md#NonlinearEigenproblems.NEPSolver.resinv),
[`contour_beyn`](methods.md#NonlinearEigenproblems.NEPSolver.contour_beyn)
and
[`sgiter`](methods.md#NonlinearEigenproblems.NEPSolver.sgiter)).
However, many NEP-algorithms
do not require a very accurate derivative.
We now show how you can make a numerical
approximation of the derivative available, if you
do not want to compute the exact derivative.
The example below uses finite differences, but
any numerical differentiation procedure may be used.
(The code does not use derivatives in `mder`,
since all calls are done with `der=0`.)

```julia
n=100;
function my_Mder_FD(λ::Float64,der::Int=0)
  if (der>1)
   error("Higher derivatives not supported");
  end
  # 3*n nnz elements in matrix
  I=Vector{Int}(undef,3*n);
  J=Vector{Int}(undef,3*n);
  F1=Vector{Float64}(undef,3*n);
  ccall(Libdl.dlsym(mylib,:mder_), Nothing,
     (Ref{Float64}, Ref{Int},Ref{Int},  Ptr{Int}, Ptr{Int}, Ptr{Float64}),
     λ, n, 0, I, J, F1)
  if (der==0)
     return sparse(I,J,F1);
  end

  if (der==1)
     # Make another fortran call to make a finite difference approximation
     ee=sqrt(eps());
     F2=Vector{Float64}(undef,3*n);
     ccall(Libdl.dlsym(mylib,:mder_), Nothing,
          (Ref{Float64}, Ref{Int},Ref{Int},  Ptr{Int}, Ptr{Int}, Ptr{Float64}),
          λ-ee, n, 0, I, J, F2)
     # We exploit the fact that the sparsity pattern is independent of λ
     Fder=(F1-F2)/ee;
     return sparse(I,J,Fder);
  end
end
```
Create the NEP and call a solver, in this case [`MSLP`](methods.md#NonlinearEigenproblems.NEPSolver.mslp).
```julia
julia> nep=Mder_NEP(n,my_Mder_FD,maxder=1);
julia> mslp(Float64,nep,λ=-1.8,logger=1);
iter 1 err:5.14547949525844e-6 λ=-1.7941228234498503
iter 2 err:6.60477516490422e-10 λ=-1.7940561772509709
iter 3 err:5.617933513637005e-16 λ=-1.794056168678654
```


## Implementation in NEP-PACK: advanced usage

The above procedure requires that sparse matrices are created
every time the NEP is accessed. This may be computationally
demanding. A common call in NEP-PACK,
is to compute the matrix vector product ``M(λ)*v``.
If the creation of the matrix ``M(λ)`` requires considerable
computation or storage, you may want to implement
the function which directly computes the matrix vector product.
This is made available to the NEP-PACK object
as follows.

Add the following to your `myproblem.f95`:

```fortran
subroutine matvec(s,n,v,x)
  real*8, intent(in) :: s
  integer*8, intent(in) :: n
  real*8, intent(in), dimension(n):: v
  real*8, intent(out), dimension(n):: x
  integer*8 :: p
  do p = 1, n
      x(p)=2*v(p)
  end do
  do p = 1, n-1
      x(p)= x(p) - v(p+1)
      x(p+1)= x(p+1) - v(p)
  end do
  x(n)=x(n)+v(1)*s*s*s;
  x(1)=x(1)-v(n)*exp(s);
end subroutine matvec
```
After recompilation of the library file
`myproblem.so`,
restarting Julia, and loading again `myproblem.so`, we
can make a matvec function available.
```julia
function my_matvec(λ,v)
   v=vec(v);  # It has to be a vector
   x=copy(v); # Allocate a vector for storage of result
   ccall(Libdl.dlsym(mylib,:matvec_), Nothing,
      (Ref{Float64}, Ref{Int}, Ptr{Float64},  Ptr{Float64}),
      λ, n, v, x)
   return x;
end
```
We can now create a `Mder_Mlincomb_NEP` which is defined from both
matrix derivative computations as well as matrix vector products (or more
generally linear combinations of derivatives).
```julia-repl
julia> nep2=Mder_Mlincomb_NEP(n,my_Mder,1,my_matvec,0);
```
The `1` and `0` specify the highest derivative available for the two functions.
We can now solve the NEP with many methods, e.g.
[`resinv`](methods.md#NonlinearEigenproblems.NEPSolver.resinv).
```julia-repl
julia> resinv(Float64,nep2,λ=-1.8,v=ones(n),logger=1);
Precomputing linsolver
iter 1 err:0.4903565024143571 λ=-1.8
iter 2 err:0.11453605256493624 λ=-1.1926857650225988
...
iter 7 err:5.834331567428063e-13 λ=-1.7940561686808254
iter 8 err:2.989922602862964e-15 λ=-1.794056168678641
```

When using NEP-solvers requiring higher derivatives,
the above procedure can also be used to compute
linear combinations of higher derivatives by implementing
a `compute_Mlincomb` which takes a matrix as input.

![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_FORTRAN1)
