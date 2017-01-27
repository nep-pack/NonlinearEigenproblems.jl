module NEPCore

    # Fundamental nonlinear eigenvalue problems
    export NEP
    # 
    export size
    export issparse
    export NoConvergenceException
    export LinSolver
    export interpolate

    # Core interfaces
    export compute_Mder
    export compute_Mlincomb
    export compute_MM

    # NEP-functions

    export compute_resnorm
    export compute_rf


    # Helper functions 
    export compute_Mlincomb_from_MM
    export compute_Mlincomb_from_Mder

    export default_errmeasure

    import Base.size  # Overload for nonlinear eigenvalue problems
    import Base.issparse  # Overload for nonlinear eigenvalue problems



    """
Determines if a method is defined for the concrete class.
False will be returned if methodname is only defined for the
abstract superclass of s.
"""
    macro method_concretely_defined(methodname,s)
        return :(length(methodswith(typeof($s), $methodname, false))>0)
    end


    ############################################
    # Default NEP functions
    #

    abstract NEP;
    

    """
    compute_Mder(nep::NEP,λ::Number,i::Integer=0)
 Computes the ith derivative of NEP evaluated in λ\\
 Usage:\\
   `compute_Mder(nep,λ)`  # Evaluate NEP in λ
"""
    function compute_Mder(nep::NEP,λ::Number,i::Integer=0)
        error("You need to provide an implementation of Mder for this NEP")
        return 0;
    end

    """
    compute_Mlincomb(nep::NEP,λ::Number,V;a=ones(size(V,2)))
 Computes the linear combination of derivatives\\
 ``Σ_i a_i M^{(i)}(λ) v_i``
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
Computes linear combination starting with derivative startder, i.e., first
column of V is multiplied by derivative startder.

The default implementation of this is slow. Overload for specific NEP
if you want efficiency (for aug_newton, IAR, ..).
"""
    function compute_Mlincomb(nep::NEP,λ,V,a::Array,startder::Integer)
        aa=[zeros(startder);a];
        VV=[zeros(size(nep,1),startder) V]; # This is typically slow since copy is needed
        return compute_Mlincomb(nep,λ,VV,a=aa)
    end
    
    """
    compute_MM(nep::NEP,S,V)
 Computes the sum ``Σ_i M_i V f_i(S)`` for a NEP,\\
 where `S` and `V` are matrices.
"""
    function compute_MM(nep::NEP,S,V)
        error("No procedure to compute MM")
    end




    ## Helper functions 
    function compute_Mlincomb_from_MM(nep::NEP,λ::Number,V,a)
        #println("Using poor-man's compute_MM -> compute_Mlincomb")
        k=size(V,2)
        S=jordan_matrix(k,λ,T=typeof(V[1])).'
        b=zeros(size(a));
        for i=1:k
            b[i]=a[i]*factorial(Float64(i-1))
        end
        W=V*diagm(b)
        z=compute_MM(nep,S,W)*eye(k,1)
        ## activate following for debugging (verify that Mder is equivalent)
        # z2=compute_Mlincomb_from_Mder(nep,λ,V,a)
        # println("should be zero:",norm(z2-z))
        return reshape(z,size(z,1))
    end

    function compute_Mlincomb_from_Mder(nep::NEP,λ::Number,V,a)
        #println("Using poor-man's compute_Mder -> compute_Mlincomb")
        z=zeros(size(nep,1))
        for i=1:size(a,1)
            z+=compute_Mder(nep,λ,i-1)*(V[:,i]*a[i])
        end
        return z
    end


    function compute_resnorm(nep::NEP,λ,v)
        return norm(compute_Mlincomb(nep,λ,reshape(v,size(nep,1),1)))
    end

"""
    compute_rf(nep::NEP,x; y=x, target=0, λ0=target,TOL=1e-12,max_iter=10)
Computes the rayleigh functional of nep, i.e., computes λ such that
   y^TM(λ)x=0.
"""
    function compute_rf(nep::NEP,x; y=x, target=0, λ0=target,
                        TOL=1e-12,max_iter=10)
        # Ten steps of scalar Newton's method
        λ=λ0;
        Δλ=Inf
        count=0
        while (abs(Δλ)>TOL) & (count<max_iter)
            count=count+1
            z1=compute_Mlincomb(nep,λ,reshape(x,nep.n,1))
            z2=compute_Mlincomb(nep,λ,reshape(x,nep.n,1),[1],1)
            Δλ=-dot(y,z1)/dot(y,z2);
            λ=λ+Δλ
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
        error("You need to provide an implementation of size for this NEP")
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
 ### NoConvergenceException
 Exeption thrown in case an iterative method does not converge\\
 `λ` = current eigenvalue approximation\\
 `v` = current eigenvector approximation\\
 `errmeasure` = The error measure of the current eigenpair approximation\\
 `msg`
"""
    type NoConvergenceException <: Exception
        "current eigenvalue approximation"
        λ
        "current eigenvector approximation"
        v
        "The error measure of the current eigenpair approximation"
        errmeasure
        msg
    end


    # To slove linear systems (This probably needs to become an abstract class)
    type LinSolver
        solve::Function
        Afact
        function LinSolver(A)
            this=new()
            this.Afact=factorize(A)
            this.solve=function foo(x;tol=eps())
                return this.Afact\x
            end
            return this;
        end   
    end

""" Returns a Jordan matrix """

    function jordan_matrix(n::Integer,λ::Number;T::Type=Float64)
        Z=λ*eye(T,n)+diagm(ones(T,n-1),1);
    end


    function default_errmeasure(nep::NEP)
        f=function (λ,v);
            compute_resnorm(nep,λ,v)
        end
        return f
    end

end  # End Module
