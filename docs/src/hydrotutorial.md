# Tutorial: Stability of parallel shear flows

## Background
Stability analysis of flows is a very important problem in fluid mechanics. Of particular interest are applications where a transition from laminar to turbulent flow is observed and this begins with an instability in the laminar flow. Considering a mean base laminar flow $$\overline{U} = \begin{pmatrix}U(y)& 0& 0\end{pmatrix}$$ and the perturbation $$u' = \begin{pmatrix}u& v& w\end{pmatrix}$$. The normal vorticity is defined as $$\eta = \dfrac{\partial u}{\partial z}-\dfrac{\partial w}{\partial x}$$. Linearizing the Navier-Stokes equations around the mean flow and then eliminating pressure gives us the Orr-Sommerfeld and Squire equations, which are a system of fourth order PDEs describing the dynamics of $v$ and $\eta$.

```math
\left(\Big(\dfrac{\partial }{\partial t}+U\dfrac{\partial }{\partial x}\Big)\nabla^2-U''\dfrac{\partial }{\partial x}-\frac{1}{Re}\nabla^4\right)v = 0,\\
\left(\dfrac{\partial }{\partial t}+U\dfrac{\partial }{\partial x}-\frac{1}{Re}\nabla^2\right)\eta = -U''\dfrac{\partial v}{\partial z}. 
```

## Formulation as a nonlinear eigenvalue problem
The ansatz $v(x,y,z,t) = \tilde{v}(y)\exp(i(\alpha x+\beta z-\omega t))$ and $\eta(x,y,z,t) = \tilde{\eta}(y)\exp(i(\alpha x+\beta z-\omega t))$, implying plane wave perturabtions transforms the system to

```math
\left((-i\omega+i\alpha U)(\mathcal{D}^2-\alpha^2-\beta^2)-i\alpha U''-\frac{1}{Re}(\mathcal{D}^2-\alpha^2-\beta^2)^2\right)\tilde{v} = 0,\\
\left((-i\omega+i\alpha U)-\frac{1}{Re}(\mathcal{D}^2-\alpha^2-\beta^2)\right)\eta = -i\beta U'\tilde{v},
```
where $\mathcal{D}$ denotes the 1-D differential operator $\frac{\partial }{\partial y}$. In this example, we consider the boundary conditions $\tilde{v} = \mathcal{D}\tilde{v} = \tilde{\eta} = 0$.

We assume that $$\beta$$ and $$\omega$$ are given and we wish to solve for $$\alpha$$. This is usually done by using a transformation of the form $$\begin{pmatrix}\tilde{v}\\ \tilde{\eta}\end{pmatrix} = \begin{pmatrix}\tilde{V}\\ \tilde{E}\end{pmatrix}\exp(-\alpha y)$$ which reduces the power of $$\alpha$$ from four to two. The problem is then discretized and solved as a quadratic eigenvalue problem. Since NEP-PACK allows us to tackle polynomial eigenvalue problems of any order, we can simply discretize $$\mathcal{D}$$ to $$D$$ and expand, which leads to a polynomial eigenvalue problem of fourth order. 

We define diagonal matrices $$U_0$$, $$U_1$$ and $$U_2$$  as follows.
```math
U_0(i,i) = U(y_i),\\
U_1(i,i) = U'(y_i),\\
U_2(i,i) = U''(y_i).\\
```
Here $$\{y_i\}_{i=1,\ldots,n}$$ denotes the y-coordinates of the $$n$$ grid points used for discretization. The discretized problem on expansion is

```math
\left[-\frac{I}{Re}\alpha^4-iU_0\alpha^3+\left(\frac{2D^2}{Re}+\left(i\omega-
\frac{2\beta^2}{Re}\right)I\right)\alpha^2+\left(iU_0(D^2-\beta^2I)-U_2\right)\alpha +\Big(\frac{2\beta^2D^2}{Re}-\frac{D^4}{Re}-\frac{\beta^4I}{Re}+i\omega(\beta^2I-D^2)\Big)\right]\tilde{v} = 0\\
i\beta U_1\tilde{v}+\left[\frac{I}{Re}\alpha^2 + iU_0 \alpha +\left(\left(\frac{\beta^2}{Re}-i\omega\right)I-\frac{D^2}{Re}\right)\right]\tilde{\eta}
```


## Implementation in NEP-PACK


