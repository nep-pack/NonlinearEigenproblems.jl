# Methods / functions in common for several contour integral methods
using Distributed, LinearAlgebra, Random


export MatrixIntegrator, MatrixTrapezoidal;

"""
    abstract type MatrixIntegrator

This type is used for integration of (matrix valued)
functions with a particular structure. It
is used in `contour_beyn` and `contour_block_SS`.

In order to specify your own way to
do numerical quadrature you should inherit from this
type
```
julia> abstract type MyIntegrator <: MatrixIntegrator; end
```
and implement the method `integrate_interval`:
```
julia> import NonlinearEigenproblems.NEPSolver.integrate_interval
julia> function integrate_interval(ST::Type{MyType},::Type{T},f,gv,a,b,N,logger) where {T<:Number}
```
The function `integrate_interval` should integrate functions on the interval `[a,b]` with `N` quadrature points. Further parameters:

* The type `T` (typically `ComplexF64`) specifies what the target number type is.
* The `logger::Logger` specifies if / how things should written to the log.
* `f::Function`, takes a scalar input (in the interval `[a,b]`) and returns a matrix or a vector.
* `gv::Vector{Function}` which takes scalar input (in the interval `[a,b]`)  and returns scalars.

The function `integrate_interval` should return a tensor `I`
where the last dimension is `m=length(gv)` and should contain the
integrals. More precisely,
 `I[:,:,1]`, ... `I[:,:,m]` should contain approximations
of the product of `f(x)gv[1](x)`,...`f(x)gv[m](x)` over
`[a,b]`, e.g., `I[:,:,j]` should contain an approximation of
```math
\\int_a^b f(x)g_j(x)\\,dx
```
Usually, `f(x)` is considerably more expensive to
evaluate than `g_1(x),..,g_m(x)`.

See also: `contour_beyn`, `contour_block_SS`, `MatrixTrapezoidal`
"""
abstract type MatrixIntegrator ; end


"""
    abstract type MatrixTrapezoidal <: MatrixIntegrator; end

Implements the trapezoidal rule to be used
with quadrature methods.

See also: `contour_beyn`, `contour_block_SS`, `MatrixIntegrator`
"""
abstract type MatrixTrapezoidal <: MatrixIntegrator; end



function integrate_interval(ST::Type{MatrixTrapezoidal},::Type{T},f,gv,a,b,N,logger) where {T<:Number}
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    f1=f(t[1]);
    m=size(gv,1);
    S = zeros(T,size(f1)...,m)

    push_info!(logger,string(ST)*": computing G",continues=true)
    # Matrix of gv evaluations
    G = zeros(T,N,m);
    for i=1:m
        push_info!(logger,".",continues=true);
        push_info!(logger,2,"t=$t",continues=true);
        gg=gv[i].(Vector(t));
        G[:,i] = gg
    end
    push_info!(logger,".",continues=false);

    # Do the sum
    push_info!(logger,string(ST)*": summing terms",continues=true);
    for i = 1:N
        push_info!(logger,".",continues=true);
        if (i==1) # Already computed
            temp = f1;
        else
            temp = f(t[i])
        end
        for j=1:m
            S[:,:,j] += temp*G[i,j];
        end
    end
    push_info!(logger,"")
    return S * h
end
