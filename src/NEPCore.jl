module NEPCore
    using SparseArrays
    using LinearAlgebra

    # Fundamental nonlinear eigenvalue problems
    export NEP
    #
    export NoConvergenceException
    export LostOrthogonalityException
    export interpolate

    # Core interfaces
    export compute_Mder
    export compute_Mlincomb
    export compute_Mlincomb!
    export compute_MM

    # NEP-functions

    export compute_resnorm
    export compute_rf


    # Helper functions  (avoid using these directly in NEP-methods)
    export compute_Mlincomb_from_MM
    export compute_Mlincomb_from_Mder
    export compute_Mder_from_MM

    export compute_Mlincomb_from_MM!

    import Base.size  # Overload for nonlinear eigenvalue problems
    import SparseArrays.issparse  # Overload for nonlinear eigenvalue problems




    ############################################
    # Default NEP functions
    #

    """
    abstract NEP
A `NEP` object represents a nonlinear eigenvalue problem. All NEPs should implement
```julia-repl
size(nep::NEP,d)
```
and at least one of the following

* M = [`compute_Mder(nep::NEP,λ::Number,i::Integer=0)`](@ref)
* V = [`compute_Mlincomb!(nep::NEP,λ::Number,V::AbstractVecOrMat,a::Vector)`](@ref)
* MM = [`compute_MM(nep::NEP,S,V)`](@ref)


"""
    abstract type NEP end


"""
    compute_Mder(nep::NEP,λ::Number [,i::Integer=0])

Computes the ith derivative of `nep` evaluated in `λ`.

# Example
This example shows that `compute_Mder(nep,λ,1)` gives the first derivative.
```julia-repl
julia> nep=nep_gallery("dep0");
julia> ϵ=1e-5;
julia> Aminus=compute_Mder(nep,λ-ϵ);
julia> Aminus=compute_Mder(nep,λ-ϵ);
julia> Aplus=compute_Mder(nep,λ+ϵ);
julia> opnorm((Aplus-Aminus)/(2ϵ)-compute_Mder(nep,λ,1))
1.990970375089371e-11
```
"""
    function compute_Mder(nep::NEP,λ::Number,i::Integer=0)
        error("You need to provide an implementation of compute_Mder for this NEP.\nIf you have a compute_MM-function you may want to define: \ncompute_Mder($(typeof(nep)),λ::Number,i::Integer)=compute_Mder_from_MM(nep,λ,i)")
    end


"""
    compute_Mlincomb(nep::NEP,λ::Number,V, a::Vector=ones(size(V,2)), startder=0)
    compute_Mlincomb!(nep::NEP,λ::Number,V, a::Vector=ones(size(V,2)), startder=0)
Computes the linear combination of derivatives\\
``Σ_i a_i M^{(i)}(λ) v_i``
starting from derivative `startder`. The function `compute_Mlincomb!`
does the same but may modify the `V` matrix/array.

# Example
This example shows that `compute_Mder` gives a result consistent with `compute_Mlincomb`. Note that `compute_Mlincomb` is in general faster since no matrix needs to be constructed.
```julia-repl
julia> nep=nep_gallery("dep0");
julia> v=ones(size(nep,1)); λ=-1+1im;
julia> norm(compute_Mder(nep,λ,1)*v-compute_Mlincomb(nep,λ,hcat(v,v),[0,1]))
1.0778315928076987e-15

```
"""
compute_Mlincomb!(nep::NEP,λ::Number,V::AbstractVecOrMat,a::Vector), compute_Mlincomb(nep::NEP,λ::Number,V::AbstractVecOrMat, a::Vector)

    function compute_Mlincomb!(nep::NEP,λ::Number,V::AbstractVecOrMat,a::Vector)
        # This will manually scale the columns in V by the vector a.
        if (ones(eltype(a),size(a,1))==a) # No scaling necessary
            return compute_Mlincomb!(nep,λ,V);
        end
        if (isa(V,AbstractVector))
            V[:]=V*a[1];
        else
            D=Diagonal(a);
            rmul!(V,D);
        end
        return compute_Mlincomb!(nep,λ,V);
    end

    # Recommend to make a copy of V and call compute_Mlincomb! if function not available
    function compute_Mlincomb(nep::NEP,λ::Number,V::AbstractVecOrMat)
        @warn "It seems you have not implemented compute_Mlincomb(nep,λ,V) for this NEPType. If you have implemented compute_Mlincomb! you need to add \ncompute_Mlincomb(nep::$(typeof(nep)),λ::Number,V::AbstractVecOrMat)=compute_Mlincomb!(nep,λ,copy(V))"
        error("No compute_Mlincomb(nep,λ,V) implemented (typeof(nep)=",typeof(nep),")")
    end
    compute_Mlincomb(nep::NEP,λ::Number,V::AbstractVecOrMat, a::Vector)=
           compute_Mlincomb!(nep,λ,copy(V), a)

    # Note: The following function is commented out since default behaviour is
    # by to manually create a bigger a-vector (and call without startder) see below
    #compute_Mlincomb(nep::NEP,λ::Number,V::AbstractVecOrMat, a::Vector, startder::Integer)=compute_Mlincomb!(nep,λ,copy(V), a, startder)


    # Default behavior of the compute_Mlincomb! is to just call compute_Mlincomb
    compute_Mlincomb!(nep::NEP,λ::Number,V::AbstractVecOrMat, a::Vector, startder::Integer)=compute_Mlincomb(nep,λ,V, a, startder)
    # Note: The following function is commented out since, default behaviour is
    # by manual scaling of columns (see above), not calling compute_Mlincomb()
    # compute_Mlincomb!(nep::NEP,λ::Number,V::AbstractVecOrMat, a::Vector)=compute_Mlincomb(nep,λ,V, a) # This is instead achieved by
    compute_Mlincomb!(nep::NEP,λ::Number,V::AbstractVecOrMat)=compute_Mlincomb(nep,λ,V)

"""
    compute_Mlincomb(nep::NEP,λ::Number,V,a::Array,startder::Integer)

Computes linear combination starting with derivative startder, i.e.,
``Σ_i a_i M^{(i+startder)}(λ) v_i``

The default implementation of this can be slow. Overload for specific NEP
if you want efficiency, e.g., in  `augnewton`, `iar`, and others.
"""
    function compute_Mlincomb(nep::NEP,λ::Number,V::AbstractVecOrMat,a::Vector,startder::Integer)
        aa=[zeros(eltype(a), startder);a];
        VV=[zeros(eltype(V), size(nep,1),startder) V]; # This is typically slow since copy is needed
        return compute_Mlincomb(nep,λ,VV,aa)
    end

"""
    compute_MM(nep::NEP,S,V)

Computes the sum ``Σ_i M_i V f_i(S)`` for a NEP, where ``S`` and ``V`` are matrices, and the NEP satisfies ``M(λ)=Σ_i M_i f_i(λ)``.

# Example
This example shows that for diagonal `S`, the result of `compute_MM` can
also be computed with `compute_Mlincomb`
```julia-repl
julia> nep=nep_gallery("dep0");
julia> D=diagm(0 => [1,2])
2×2 Array{Int64,2}:
 1  0
 0  2
julia> V=ones(size(n,1),2);
julia> W=compute_MM(nep,D,V);
julia> norm(W[:,1]-compute_Mlincomb(nep,D[1,1],V[:,1]))
1.1102230246251565e-16
julia> norm(W[:,2]-compute_Mlincomb(nep,D[2,2],V[:,2]))
0.0
```
# Reference
Properties of the quantity ``Σ_i M_i V f_i(S)`` for
non-polynomial nonlinear eigenvalue problems were
extensively used in:
* D. Kressner A block Newton method for nonlinear eigenvalue problems, Numer. Math., 114 (2) (2009), pp. 355-372
* C. Effenberger, Robust solution methods for nonlinear eigenvalue problems, PhD thesis, 2013, EPF Lausanne

"""
    function compute_MM(nep::NEP,S,V)
        error("No procedure to compute MM (typeof(nep)=",typeof(nep),")")
    end




    ## Helper functions
    """
    compute_Mlincomb_from_MM(nep::NEP,λ::Number,V,a)
The function computes Mlincomb by a call to compute_MM. The relationship between Mlincomb and MM is described in issue #2 and #3
Usage normally by overloading:
    compute_Mlincomb(nep::MyNEP,λ::Number,V,a)=compute_Mlincomb_from_MM(nep,λ,V,a)
"""
    compute_Mlincomb_from_MM(nep::NEP,λ::Number,V,a)=compute_Mlincomb_from_MM!(nep,λ,copy(V),copy(a))
    """
    compute_Mlincomb_from_MM!(nep::NEP,λ::Number,V,a)
Same as [`compute_Mlincomb`](@ref), but modifies V and a.
"""
    function compute_Mlincomb_from_MM!(nep::NEP,λ::Number,V,a::Array{<:Number,1})
        # This function it is based on the evaluation of matrix function of a bidiagonal matrix
        # Should we document the methematical
        k=size(V,2);
        # we need to assume that the elements of a are different than zero.
        V[:,findall(x->x==0,a)] .= 0
        a[findall(x->x==0,a)] .= 1
        S=diagm(0 => λ*ones(eltype(V),k)) + diagm(-1 => (a[2:k]./a[1:k-1]).*(1:k-1))
        z=compute_MM(nep,S,V)[:,1];
        return a[1]*reshape(z,size(z,1))
    end
    """
    compute_Mlincomb_from_Mder(nep::NEP,λ::Number,V,a)
The function computes Mlincomb by a call to compute_Mder. This function is slow since it requires the construction of the matrices.
Usage normally by overloading:
    compute_Mlincomb(nep::MyNEP,λ::Number,V,a)=compute_Mlincomb_from_Mder(nep,λ,V,a)
"""
    function compute_Mlincomb_from_Mder(nep::NEP,λ::Number,V,a::Array{<:Number,1})
        #println("Using poor-man's compute_Mder -> compute_Mlincomb")
        z=zeros(size(nep,1))
        for i=1:length(a)
            if (a[i] != 0)
                z+=compute_Mder(nep,λ,i-1)*(V[:,i]*a[i])
            end
        end
        return z
    end

    """
    compute_Mder_from_MM(nep::NEP,λ::Number,i::Integer=0)
Computes the Mder function from MM using the fact that MM of
a jordan block becomes derivatives
"""
    function compute_Mder_from_MM(nep::NEP,λ::Number,i::Integer=0)
        J=transpose(jordan_matrix(typeof(λ),i+1,λ))
        n=size(nep,1);
        S=kron(J, Matrix(1.0I, n, n))
        V=factorial(i) * kron(sparse(1.0I, 1, i+1)[:,end:-1:1], sparse(1.0I, n, n))
        W=compute_MM(nep,S,V)
        return W[1:n,1:n]
    end

    """
    compute_resnorm(nep::NEP,λ,v)
Computes the residual norm of the `nep`, in the point `λ`, with the vector
`v`, i.e., ``||M(λ)v||``.
"""
    function compute_resnorm(nep::NEP,λ,v)
        return norm(compute_Mlincomb(nep,λ,reshape(v,size(nep,1),1)))
    end

  """
    compute_rf([eltype],nep::NEP,x; y=x, target=zero(T), λ0=target,TOL=eps(real(T))*1e3,max_iter=10)

Computes the Rayleigh functional of nep, i.e., computes a vector ``Λ`` of values ``λ``
such that ``y^TM(λ)x=0``. The default behaviour consists of a scalar valued
Newton-iteration, and the returned vector has only one element.

The given eltype<:Number is the type of the returned vector.

# Example

```julia-repl
julia> nep=nep_gallery("dep0");
julia> x=ones(size(nep,1));
julia> s=compute_rf(Float64,nep,x)[1]; # Take just first element
0.6812131933795569
julia> x'*compute_Mlincomb(nep,s,x)
-8.881784197001252e-16
```
"""
    compute_rf(nep::NEP,x;params...) = compute_rf(ComplexF64,nep,x;params...)
    function compute_rf(::Type{T}, nep::NEP, x; y=x, target=zero(T), λ0=target,
                        TOL=eps(real(T))*1e3, max_iter=10) where T
        # Ten steps of scalar Newton's method
        λ_iter = T(λ0);
        Δλ = T(Inf)
        count = 0
        while (abs(Δλ)>TOL) & (count<max_iter)
            count = count+1
            z1 = compute_Mlincomb(nep, λ_iter, reshape(x,size(nep,1),1))
            z2 = compute_Mlincomb(nep, λ_iter, reshape(x,size(nep,1),1),[T(1)],1)

            Δλ =- dot(y,z1)/dot(y,z2);
            λ_iter += Δλ
        end

        # Return type is a vector of correct type
        λ_star::Array{T,1} = Array{T,1}(undef, 1)
        if (T <: Real) && (typeof(λ_iter) != T) && (imag(λ_iter)/real(λ_iter) < TOL)
            # Looking for a real quantity (AND) iterate is not real (AND) complex part is negligible
            λ_star[1] = real(λ_iter) # Truncate to real
        else
            λ_star[1] = λ_iter
        end
        return λ_star
    end


   """
    size(nep::NEP)
    size(nep::NEP,dim)
Overloads the size functions for NEP.\\
Size returns the size of the matrix defining the NEP.

Note: All NEPs must implement this function.
"""
    function size(nep::NEP)
        error("You need to provide an implementation of size for this NEP.")
    end
    function size(nep::NEP,dim)
        error("You need to provide an implementation of size for this NEP.")
    end


   """
    issparse(nep::NEP)
Overloads the issparse functions for NEP.\\
Issparse returns `true` if the undelying type of the NEP is\\
sparse, and `false` if it is dense.\\
Default behaviour: Check sparsity of `compute_Mder(nep,0)`

"""
    function issparse(nep::NEP)
        issparse(compute_Mder(nep,0.0)) # TODO: This might require a redesign when/if NEPs are parametric. Type of 0.0?
    end





    ############################################
    # Misc helpers
    #

    """
    struct NoConvergenceException
Exeption thrown in case an iterative method does not converge\\
`λ` = current eigenvalue(s) approximation\\
`v` = current eigenvector(s) approximation\\
`errmeasure` = The error measure of the current eigenpair(s) approximation\\
`msg`
"""
    struct NoConvergenceException <: Exception
        "current eigenvalue(s) approximation"
        λ
        "current eigenvector(s) approximation"
        v
        "The error measure of the current eigenpair(s) approximation"
        errmeasure
        msg
    end
    # Avoid dumping the huge eigenvectors in case of exception
    Base.showerror(io::IO, e::NoConvergenceException) =
        print(io, "No convergence: '",e.msg,"' ",
              "eigenvalue approx:",e.λ,", errmeasure:",e.errmeasure)


""" Returns a Jordan matrix """
    jordan_matrix(n::Integer,λ::Number)=jordan_matrix(ComplexF64,n,λ)
    function jordan_matrix(::Type{T},n::Integer,λ::Number) where T<:Number
        Z = T(λ) * Matrix{T}(I, n, n) + diagm(1 => ones(T, n-1))
    end


    """
    struct LostOrthogonalityException
`msg`
"""
    struct LostOrthogonalityException <: Exception
        msg
    end


end  # End Module
