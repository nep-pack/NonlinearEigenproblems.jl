# Tutorial: Application to absorbing boundary conditions
## A Schrödinger equation

Consider Schrödinger type eigenvalue problem on the interval $[0,L_1]$,
```math
\begin{eqnarray*}
 \left(
 \frac{\partial^2}{\partial x^2}
-V(x)-\lambda^2
\right)\psi(x)&=&0,\,\; x\in[0,L_1]\\
   \psi(0)&=&0\\
   \psi(L_1)&=&0.
\end{eqnarray*}
```
We wish to compute eigenvalue ``λ`` and eigenfunction ``\psi``.
Moreover, we assume that the potential function ``V(x)`` is benign in
the domain ``[L_0,L_1]``, e.g., constant. In the simulations we will consider this
function
```math
  V(x)=
\begin{cases}
1+\sin(\alpha x)  & x\in[0,L_0]=[0,1]\\
V_0 & x\in(L_0,L_1)=(1,8)
\end{cases}
```
where ``α`` is large, i.e., the potential has high oscillations
in one part of the domain.
This tutorial illustrates how we can avoid a discretization
of the domain ``[L_0,L_1]`` and only discretize ``[0,L_0]``,
by solving a NEP.

## Derivation of reduced domain differential equation

The technique is based on moving the boundary condition at ``L_1``
to ``L_0``. This can be done without doing any approximation,
if we allow the new artificial boundary condition at ``L_1``
to depend on ``λ``. The technique is called absorbing boundary conditions.



We first note that we can transform the problem to first order form
```math
  \frac{d}{dx}
\begin{bmatrix}\psi(x)\\\psi'(x)\end{bmatrix}
=
\begin{bmatrix}
0 & 1\\
\lambda^2+V(x) & 0
\end{bmatrix}
\begin{bmatrix}\psi(x)\\\psi'(x)\end{bmatrix},
```
We now apply this in the domain ``[L_0,L_1]``,
the fact that ``V(x)`` is constant allows us to directly
express the solution using the matrix exponential
```math
\begin{bmatrix}\psi(x)\\\psi'(x)\end{bmatrix}
=\exp\left((x-L_0)
\begin{bmatrix}
0 & 1\\
\lambda^2+V_0 & 0
\end{bmatrix}
\right)
\begin{bmatrix}\psi(L_0)\\\psi'(L_0)\end{bmatrix}.
```
when ``x\in[L_0,L_1]``. The boundary condition ``\psi(L_1)=0`` can be imposed
as
```math
0=\psi(L_1)=\begin{bmatrix}1 & 0\end{bmatrix}
\begin{bmatrix}\psi(L_1)\\\psi'(L_1)\end{bmatrix}
=\begin{bmatrix}1 & 0\end{bmatrix}\exp\left((L_1-L_0)
\begin{bmatrix}
0 & 1\\
\lambda^2+V_0 & 0
\end{bmatrix}
\right)
\begin{bmatrix}\psi(L_0)\\\psi'(L_0)\end{bmatrix}.
```
By explicitly computing the matrix exponential of the antidiagonal two-by-two
matrix we obtain the relation
```math
0=
g(λ)\psi(L_0)+
f(λ)\psi'(L_0).
```
where
```math
g(λ):=\sqrt{λ^2+V_0}\cosh\left((L_1-L_0)\sqrt{λ^2+V_0}\right)
```
```math
f(λ):=\sinh\left((L_1-L_0)\sqrt{λ^2+V_0}\right).
```
We now observe that the differential equation
for the domain ``[0,L_0]`` is
```math
\begin{eqnarray*}
 \left(
 \frac{\partial^2}{\partial x^2}
-V(x)-\lambda^2
\right)\psi(x)&=&0,\,\; x\in[0,L_0]\\
   \psi(0)&=&0\\
   g(λ)\psi(L_0)+f(λ)\psi'(L_1)&=&0.
\end{eqnarray*}
```
which is a boundary value problem
with a mixed boundary condition (since it contains
both ``\psi(L_0)`` and ``\psi'(L_0)``).


## Discretization of the λ-dependent boundary value problem

The λ-dependent boundary value problem can now be discretized in the domain ``[0,L_0]``
with usual techniques, e.g., with finite difference as follows.

Let $x_k=hk$, $k=1,\ldots n$ and $h=1/n$ such that $x_1=h$ and
$x_n=1=L_0$. Approximation of the $\lambda$-dependent boundary condition
can be achieved with the a one-sided second order
difference scheme
```math
   0=g(λ)\psi(L_0)+f(λ)\frac{1}{h}\left(\frac32 \psi(L_0)
-2\psi(x_{n-1})
+\frac12\psi(x_{n-2})\right)+O(h^2)
```
Let
```math
  D_n=
\frac1{h^2}
\begin{bmatrix}
-2  & 1 & 0 &\\
1 & \ddots &1& \\
0 & 1 &-2 & 1\\
0 & \cdots & 0 & 0
\end{bmatrix}\textrm{ and }
\underline{I}_n=\begin{bmatrix}1 &\\ &\ddots\\ &&1 \\ 0&\cdots &0\end{bmatrix}
```
Then the boundary value problem can expressed as ``M(λ)v=0`` where
```math
M(λ)=D_n-\operatorname{diag}(V(x_1),\ldots,V(x_{n-1}),0)-λ^2\underline{I}_n
+g(λ)e_ne_n^T+f(λ)F
```
and
```math
F=\frac{1}{2h}e_ne_{n-2}^T-\frac{2}{h}e_ne_{n-1}^T+\frac{3}{2h}e_ne_n^T
```

## Implementation in NEP-PACK
The above discretization can be expressed as a [`SPMF`](types.md#SPMF) with
four terms. Let us set up the matrices first
```julia
L0=1; L1=8; V0=10.0;
xv=Vector(range(0,stop=L0,length=1000))
h=xv[2]-xv[1];
n=size(xv,1);
α=25*pi/2;
V=x->1+sin(α*x);
Dn=spdiagm(-1 => [ones(n-2);0]/h^2, 0 => -2*ones(n-1)/h^2, 1 => ones(n-1)/h^2)
Vn=spdiagm(0 => [V.(xv[1:end-1]);0]);
In=spdiagm(0 => [ones(n-1);0])
F=sparse([n, n, n],[n-2, n-1, n],[1/(2*h), -2/h, 3/(2*h)])
G=sparse([n],[n],[1]);
```
The corresponding functions
```julia
f1=S->one(S);
f2=S->-S^2;
hh=S-> sqrt(S^2+V0*one(S))
g=S-> cosh((L1-L0)*hh(S))
f=S-> inv(hh(S))*sinh((L1-L0)*hh(S))
```
Note that `sinh(A)` and `cosh(A)` for matrices `A`, are interpreted as matrix functions (not elementwise) which is also what we need. Create the NEP and solve it:
```julia
using NonlinearEigenproblems
nep=SPMF_NEP([Dn-Vn,In,G,F],[f1,f2,g,f]);
```

We can apply a NEP-solver and compute solutions
```julia
(λ1,v1)=quasinewton(nep,displaylevel=1,λ=2.5im,v=ones(n),tol=1e-9);
(λ2,v2)=quasinewton(nep,displaylevel=1,λ=3.2im,v=ones(n),tol=1e-9)
(λ3,v3)=quasinewton(nep,displaylevel=1,λ=5im,v=ones(n),tol=1e-9)
(λ4,v4)=quasinewton(nep,displaylevel=1,λ=6im,v=ones(n),tol=1e-9)
```
We can plot the solutions:
```julia
using Plots
plot(xv,real(v1)/norm(v1))
plot!(xv,real(v2)/norm(v2))
plot!(xv,real(v3)/norm(v3))
plot!(xv,real(v4)/norm(v4))
```
resulting in
