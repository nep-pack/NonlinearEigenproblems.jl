# Tutorial: Hydrodynamic stability as a fourth order PEP

## Background
Stability analysis of flows is a very important problem in fluid mechanics. Of particular interest are applications where a transition from laminar to turbulent flow is observed and this begins with an instability in the laminar flow. Considering a mean base laminar flow $$\overline{U} = \begin{pmatrix}U(y)& 0& 0\end{pmatrix}$$ and the perturbation $$u' = \begin{pmatrix}u& v& w\end{pmatrix}$$. Linearizing the Navier-Stokes equations around the mean flow gives us the Orr-Sommerfeld and Squire equations

```math
\dfrac{\partial u}{\partial t}+U\dfrac{\partial u}{\partial x}+vU' = -\dfrac{\partial p}{\partial x}+\frac{1}{Re}\nabla^2u\\
\dfrac{\partial v}{\partial t}+U\dfrac{\partial v}{\partial x} = -\dfrac{\partial p}{\partial y}+\frac{1}{Re}\nabla^2v\\
\dfrac{\partial w}{\partial t}+U\dfrac{\partial w}{\partial x} = -\dfrac{\partial p}{\partial z}+\frac{1}{Re}\nabla^2w
```
Taking 