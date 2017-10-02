module NEPCore
    # Fundamental nonlinear eigenvalue problems
    export NEP
    #
    export size
    export issparse
    export NoConvergenceException
    export LostOrthogonalityException
    export LinEigSolver
    export interpolate

    # Core interfaces
    export compute_Mder
    export compute_Mlincomb
    export compute_MM

    # NEP-functions

    export compute_resnorm
    export compute_rf


    # Helper functions  (avoid using these directly in NEP-methods)
    export compute_Mlincomb_from_MM
    export compute_Mlincomb_from_Mder
    export compute_Mder_from_MM

    export compute_Mlincomb_from_MM!

    export default_errmeasure

    import Base.size  # Overload for nonlinear eigenvalue problems
    import Base.issparse  # Overload for nonlinear eigenvalue problems



    """
Determines if a method is defined for the concrete class.
False will be returned if methodname is only defined for the
abstract superclass of s.
"""
    macro method_concretely_defined(methodname,s)
        return :( length(methodswith(typeof($(esc(s))), $(esc(methodname)), false))>0 )
    end
#    macro method_concretely_defined(methodname,s)
#        return :(length(methodswith(typeof($s), $methodname, false))>0)
#    end


    ############################################
    # Default NEP functions
    #

    """
    abstract NEP
NEP represents a nonlinear eigenvalue problem
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
julia> norm((Aplus-Aminus)/(2ϵ)-compute_Mder(nep,λ,1))
1.990970375089371e-11
```
"""
    function compute_Mder(nep::NEP,λ::Number,i::Integer=0)

        extra_msg=""
        if (@method_concretely_defined(compute_MM,nep))
            extra_msg=", or choose to call the compute_Mder_from_MM, which can be slow"
        end

        error("You need to provide an implementation of compute_Mder for this NEP, or choose to use compute_Mder_from_MM"*extra_msg*".")
        return 0;
    end

    """
    compute_Mlincomb(nep::NEP,λ::Number,V;a=ones(size(V,2)))
Computes the linear combination of derivatives\\
``Σ_i a_i M^{(i)}(λ) v_i``

# Example
This example shows that `compute_Mder` gives a result consistent with `compute_Mlincomb`. Note that `compute_Mlincomb` is in general faster since no matrix needs to be constructed.
```julia-repl
julia> nep=nep_gallery("dep0");
julia> v=ones(size(nep,1)); λ=-1+1im;
julia> norm(compute_Mder(nep,λ,1)*v-compute_Mlincomb(nep,λ,hcat(v,v),a=[0,1]))
1.0778315928076987e-15

```
"""
    function compute_Mlincomb(nep::NEP,λ::Number,V;a=ones(size(V,2)))

        # determine a default behavior (may lead to loss of performance)
        if (@method_concretely_defined(compute_MM,nep))
            return compute_Mlincomb_from_MM(nep,λ,V,a)
        elseif (@method_concretely_defined(compute_Mder,nep))
            return compute_Mlincomb_from_Mder(nep,λ,V,a)
        else
            error("No procedure to compute Mlincomb")
        end
    end

"""
    compute_Mlincomb(nep::NEP,λ::Number,V,a::Array,startder::Integer)
Computes linear combination starting with derivative startder, i.e.,
``Σ_i a_i M^{(i+startder)}(λ) v_i``

The default implementation of this can be slow. Overload for specific NEP
if you want efficiency (for aug_newton, IAR, ..).
"""
    function compute_Mlincomb(nep::NEP,λ,V,a::Array{<:Number,1},startder::Integer)
        aa=[zeros(startder);a];
        VV=[zeros(size(nep,1),startder) V]; # This is typically slow since copy is needed
        return compute_Mlincomb(nep,λ,VV,a=aa)
    end

"""
    compute_MM(nep::NEP,S,V)
Computes the sum ``Σ_i M_i V f_i(S)`` for a NEP, where `S` and `V` are matrices.

# Example
This example shows that for diagonal `S`, the result of `compute_MM` can
also be computed with `compute_Mlincomb`
```julia-repl
julia> nep=nep_gallery("dep0");
julia> D=diagm([1,2])
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
"""
    function compute_MM(nep::NEP,S,V)
        error("No procedure to compute MM")
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
        # we need to assume that the elements of a are is different than zero.
        V[:,find(x->x==0,a)]=0; a[find(x->x==0,a)]=1;
        S=diagm(λ*ones(eltype(V),k))+diagm((a[2:k]./a[1:k-1]).*(1:k-1),1); S=S.';
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
        J=sparse(jordan_matrix(typeof(λ),i+1,λ).')
        n=size(nep,1);
        S=kron(J,speye(n))
        V=factorial(i)*kron(speye(1,i+1)[:,end:-1:1],speye(n))
        W=compute_MM(nep,S,V)
        return W[1:n,1:n]
    end
    """
    compute_resnorm(nep::NEP,λ,v)
Computes the residual norm ||M(λ)v|| of the nep.
"""
    function compute_resnorm(nep::NEP,λ,v)
        return norm(compute_Mlincomb(nep,λ,reshape(v,size(nep,1),1)))
    end

  """
    compute_rf{T}(::Type{T},nep::NEP,x; y=x, target=zero(T), λ0=target,
                        TOL=eps(real(T))*1e3,max_iter=10)
Computes the rayleigh functional of nep, i.e., computes λ such that
   y^TM(λ)x=0.
"""
    compute_rf(nep::NEP,x;params...) = compute_rf(Complex128,nep,x;params...)
    function compute_rf{T}(::Type{T},nep::NEP,x; y=x, target=zero(T), λ0=target,
                        TOL=eps(real(T))*1e3,max_iter=10)
        # Ten steps of scalar Newton's method
        λ = T(λ0);
        Δλ = T(Inf)
        count = 0
        while (abs(Δλ)>TOL) & (count<max_iter)
            count=count+1
            z1=compute_Mlincomb(nep,λ,reshape(x,size(nep,1),1))
            z2=compute_Mlincomb(nep,λ,reshape(x,size(nep,1),1),[1],1)
            println("typeof x:",typeof(x)," z1:",typeof(z1)," z2:",typeof(z2),
                    " λ:",λ)
            Δλ=- T(dot(y,z1)/dot(y,z2));
            λ += Δλ
        end
        return λ
    end


   """
    size(nep::NEP,dim=-1)
Overloads the size functions for NEP.\\
Size returns the size of the matrix defining the NEP.
Note: All NEPs must implement this function.
"""
    function size(nep::NEP,dim=-1)
        error("You need to provide an implementation of size for this NEP.")
        return 0;
    end


   """
    issparse(nep::NEP)
Overloads the issparse functions for NEP.\\
Issparse returns `true` if the undelying type of the NEP is\\
sparse, and `false` if it is dense.\\
Default behaviour: Check sparsity of `compute_Mder(nep,0)`

"""
    function issparse(nep::NEP)
        issparse(compute_Mder(nep,0))
    end





    ############################################
    # Misc helpers
    #

    """
    type NoConvergenceException
Exeption thrown in case an iterative method does not converge\\
`λ` = current eigenvalue(s) approximation\\
`v` = current eigenvector(s) approximation\\
`errmeasure` = The error measure of the current eigenpair(s) approximation\\
`msg`
"""
    type NoConvergenceException <: Exception
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
    jordan_matrix(n::Integer,λ::Number)=jordan_matrix(Complex128,n,λ)
    function jordan_matrix{T<:Number}(::Type{T},n::Integer,λ::Number)
        Z=T(λ)*eye(T,n)+diagm(ones(T,n-1),1);
    end

    """
    default_errmeasure(nep::NEP)
The default way of measuring error (residual norm).
"""
    function default_errmeasure(nep::NEP)
        f=function (λ,v);
            compute_resnorm(nep,λ,v)/norm(v)
        end
        return f
    end

    """
    type LostOrthogonalityException
`msg`
"""
    type LostOrthogonalityException <: Exception
        msg
    end

end  # End Module
