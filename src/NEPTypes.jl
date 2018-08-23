
module NEPTypes
    # Specializalized NEPs
    export ProjectableNEP
    export DEP
    export PEP
    export REP
    export SPMF_NEP
    export AbstractSPMF
    export SumNEP, SPMFSumNEP, GenericSumNEP

    export Proj_NEP;
    export Proj_SPMF_NEP;
    export create_proj_NEP;

    export interpolate
    export interpolate_cheb

    export set_projectmatrices!;

    using ..NEPCore
    using PolynomialZeros
    using Polynomials

    # We overload these
    import ..NEPCore.compute_Mder
    import ..NEPCore.compute_Mlincomb
    import ..NEPCore.compute_Mlincomb!
    import ..NEPCore.compute_MM
    import ..NEPCore.compute_resnorm
    import ..NEPCore.compute_rf

    import Base.size
    import Base.issparse


    export compute_Mder
    export compute_Mlincomb
    export compute_Mlincomb!
    export compute_MM
    export compute_resnorm
    export compute_rf
    export size
    export companion
    export issparse

    export get_Av
    export get_fv




    #
    """
    abstract ProjectableNEP <: NEP

A ProjectableNEP is a NEP which can be projected, i.e., one can construct the problem W'*M(λ)Vw=0 with the Proj_NEP. A NEP which is of this must have the function `create_proj_NEP(orgnep::ProjectableNEP)` implemented. This function must return a `Proj_NEP` See also `set_projectmatrices!".

# Example:
julia> nep=nep_gallery("dep0");
julia> typeof(nep)
NEPTypes.DEP
julia> isa(nep,ProjectableNEP)
true
julia> projnep=create_proj_NEP(nep)
NEPTypes.Proj_SPMF_NEP(NEPTypes.DEP(5, Array{Float64,2}[[0.679107 0.297336 … -0.187573 -0.117138; 0.828413 0.0649475 … -1.60726 -0.601254; … ; -0.134854 -0.51421 … 2.27623 -0.0886163; 0.586617 1.57433 … 0.219693 0.279466], [0.111422 1.42305 … 0.481556 -0.185424; -0.357884 0.408387 … -0.321943 1.26972; … ; 0.300234 -0.296278 … -0.178789 -0.0671867; -0.762677 0.691111 … -1.47788 0.577282]], [0.0, 1.0]), #undef, #undef, #undef, Array{Float64,2}[[1.0 0.0 … 0.0 0.0; 0.0 1.0 … 0.0 0.0; … ; 0.0 0.0 … 1.0 0.0; 0.0 0.0 … 0.0 1.0], [0.679107 0.297336 … -0.187573 -0.117138; 0.828413 0.0649475 … -1.60726 -0.601254; … ; -0.134854 -0.51421 … 2.27623 -0.0886163; 0.586617 1.57433 … 0.219693 0.279466], [0.111422 1.42305 … 0.481556 -0.185424; -0.357884 0.408387 … -0.321943 1.26972; … ; 0.300234 -0.296278 … -0.178789 -0.0671867; -0.762677 0.691111 … -1.47788 0.577282]], Function[NEPTypes.#1, NEPTypes.#2, NEPTypes.#3])
julia> set_projectmatrices!(projnep,eye(5,1),eye(5,1))
NEPTypes.SPMF_NEP(1, Array{Float64,2}[[1.0], [0.679107], [0.111422]], Function[NEPTypes.#1, NEPTypes.#2, NEPTypes.#3], false, 1×1 SparseMatrixCSC{Float64,Int64} with 0 stored entries)

"""
    abstract type ProjectableNEP <: NEP end

    #######################################################
    ### Sum of products matrices and functions

    """
    abstract  AbstractSPMF <: ProjectableNEP

An AbstractSPMF is an abstract class representing NEPs which can be represented
as a Sum of products of matrices and functions ``M(λ)=Σ_i A_i f_i(λ)``,
where i = 0,1,2,..., all of the matrices are of size n times n and f_i are functions.
Any AbstractSPMF has to have implementations of get_Av() and get_fv() which return the
functions and matrices.
"""
    abstract  type AbstractSPMF <: ProjectableNEP end # See issue #17

    """
    get_Av(nep::AbstractSPMF)
Returns an array of matrices A_i in the AbstractSPMF: ``M(λ)=Σ_i A_i f_i(λ)``
"""
    function get_Av(nep::AbstractSPMF) # Dummy function which enforces that you have to implement
        error("You need to implement get_Av for all AbstractSPMFs")
    end
    """
    get_Av(nep::AbstractSPMF)
Returns an Array of functions (matrix functions) f_i in the AbstractSPMF: ``M(λ)=Σ_i A_i f_i(λ)``
"""
    function get_fv(nep::AbstractSPMF) # Dummy function which enforces that you have to implement
        error("You need to implement get_fv for all AbstractSPMFs")
    end

    """
    struct SPMF_NEP <: AbstractSPMF

An SPMF_NEP is a NEP defined by a Sum of Products of Matrices and Functions,
i.e.,
```math
M(λ)=∑_i A_i f_i(λ).
```
All of the matrices ``A_0,...`` are of size ``n×n``
and ``f_i`` are a functions. The  functions ``f_i`` must be defined
for matrices in the standard matrix function sense.
"""
    struct SPMF_NEP <: AbstractSPMF
         n::Int
         A::Array   # Array of Array of matrices
         fi::Array  # Array of functions
         Schur_factorize_before::Bool # Tells if you want to do the Schur-factorization at the top-level of calls to compute_MM(...)

         # Sparse zero matrix to be used for sparse matrix creation
         Zero::SparseMatrixCSC
    end
"""
     SPMF_NEP(AA,fii,Schur_fact=false)

Creates a SPMF_NEP consisting of matrices `AA` and functions `fii`. `fii` must
be an array of functions defined for matrices. `AA` is an array of
matrices. `Schur_fact` specifies if the computation of `compute_MM` should be
done by first pre-computing a Schur-factorization (which can be faster).


```julia-repl
julia> A0=[1 3; 4 5]; A1=[3 4; 5 6];
julia> id_op=S -> eye(S)
julia> exp_op=S -> expm(S)
julia> nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
julia> compute_Mder(nep,1)-(A0+A1*exp(1))
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0
```
"""
     function SPMF_NEP(AA::Vector{<:AbstractMatrix}, fii::Vector{<:Function}, Schur_fact = false)

         if (size(AA,1)==0)
             return SPMF_NEP(0); # Create empty SPMF_NEP.
         end
         n=size(AA[1],1);

         if     (size(AA,1) != 1) && (size(AA,2) == 1) # Stored as column vector - do nothing
         elseif (size(AA,1) == 1) && (size(AA,2) == 1) # It is a single entry - do nothing
         elseif (size(AA,1) == 1) && (size(AA,2) != 1) # Stored as a row-vector
             AA = vec(AA)
         else
             error("The given array should be a vector but is of size ", size(AA), ".")
         end


         if(size(AA,1) != size(fii,1))
             error("Inconsistency: Number of supplied matrices = ", length(AA), " but the number of supplied functions are = ", length(fii))
         end
         for i = 2:size(AA,1)
             if size(AA[i]) != size(AA[1])
                 error("The dimensions of the matrices mismatch: size(AA[1])",
                       size(AA[1]),"!=",size(AA[i]),"=size(AA[",i,"])")
             end
         end


         if (issparse(AA[1]))
             Zero=spones(AA[1]);
             for i=2:size(AA,1)
                 Zero=Zero+spones(AA[i]);
             end
             Zero=(Zero*1im)*0
         else
             Zero=zeros(n,n)
         end


         this=SPMF_NEP(n,AA,fii,Schur_fact,Zero);
         return this
    end
    function SPMF_NEP(n) # Create an empty NEP of size n x n
         Z=zeros(n,n)
         return SPMF_NEP(n,Vector{Matrix}(),Vector{Function}(),false,Z);
    end
    function compute_MM(nep::SPMF_NEP,S,V)
        if (issparse(V))
            if (size(V)==size(nep))
                # Initialize with zero sparse matrix which
                # has sparsity pattern already consistent
                # with sparsity pattern of M() for optimization
                Z=copy(nep.Zero)
            else
                Z=spzeros(eltype(V),size(V,1),size(V,2))
            end
        else
            Z=zeros(eltype(V),size(V))
        end
        # Sum together all the terms in the SPMF:
        if(nep.Schur_factorize_before) #Optimize if allowed to factorize before
            (T, Q, ) = schur(S)
        end
        for i=1:size(nep.A,1)
            ## Compute Fi=f_i(S) in an optimized way
            if (isdiag(S)) # optimize if S is diagonal
                Sd=diag(S);
                if (norm(Sd-Sd[1])==0) # Optimize further if S is a
                                       # multiple of identity
                    Fid=nep.fi[i](reshape([Sd[1]],1,1))[1]*ones(size(Sd,1))
                else  # Diagonal but not constant
                    Fid=zeros(Complex128,size(S,1))
                    for j=1:size(S,1)
                        Fid[j]=nep.fi[i](reshape([Sd[j]],1,1))[1]
                    end
                end
                Fi=spdiagm(Fid);
            else  # Otherwise just compute the matrix function operation
                if(nep.Schur_factorize_before)
                    Fi= Q*nep.fi[i](T)*Q'
                else
                    Fi=nep.fi[i](S)
                end
            end
            ## Sum it up
            Z=Z+nep.A[i]*(V*Fi);
        end
        return Z
    end
    function compute_Mder(nep::SPMF_NEP,λ::Number,i::Integer=0)
        if (i==0)
            Z=copy(nep.Zero)
            for i=1:size(nep.A,1)
                Z=Z+nep.A[i]*nep.fi[i](reshape([λ],1,1))[1]
            end
            return Z
        else
            # This is typically slow for i>1 (can be optimized by
            # treating the SPMF-terms individually)
            return compute_Mder_from_MM(nep,λ,i)
        end
    end
    # Use MM to compute Mlincomb for SPMFs
    compute_Mlincomb(nep::SPMF_NEP,λ::Number,
                     V::Union{AbstractMatrix,AbstractVector},a::Vector=ones(size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)

    ###########################################################
    # Delay eigenvalue problems - DEP
    #

    """
    type DEP <: AbstractSPMF
Delay eigenvalue problem
A DEP (Delay Eigenvalue problem) is defined by
the sum  ``-λI + Σ_i A_i exp(-tau_i λ)`` where all
of the matrices are of size n times n.\\
Constructor: DEP(AA,tauv) where AA is an array of the
```math
\\frac{n!}{k!(n - k)!} = \\binom{n}{k}
```
matrices A_i, and tauv is a vector of the values tau_i
"""
    type DEP{T<:AbstractMatrix} <: AbstractSPMF
        n::Int
        A::Array{T,1}     # An array of matrices (full or sparse matrices)
        tauv::Array{Float64,1} # the delays
    end
    function DEP(AA::Array{T,1},tauv::Vector=[0,1.0]) where {T<:AbstractMatrix}
        n=size(AA[1],1)
        tauvconv=Vector{Float64}(tauv);
        if (real(eltype(AA[1])) != Float64)
            warning("The delay in a DEP has hardcoded type Float64");
        end

        this=DEP{T}(n,AA,tauvconv);
        return this;
    end

# Compute the ith derivative of a DEP
    function compute_Mder(nep::DEP,λ::Number,i::Integer=0)
        local M,I;
        # T is eltype(nep.A[1]) unless λ complex, then T is complex(eltype(nep.A[1]))
        # T can be determined compile time, since DEP parametric type
        T=isa(λ,Complex)?complex(eltype(nep.A[1])) : eltype(nep.A[1]);

        if (issparse(nep.A[1])) # Can be determined compiled time since DEP parametric type
            M=spzeros(T,nep.n,nep.n)
            I=speye(T,nep.n,nep.n)
        else
            M=zeros(T,nep.n,nep.n)
            I=eye(T,nep.n,nep.n)
        end
        if i==0; M=-λ*I;  end
        if i==1; M=-I; end
        for j=1:size(nep.A,1)
            a=exp(-nep.tauv[j]*λ)*(-nep.tauv[j])^i;
            M += nep.A[j]*a
        end
        return M
    end




#  Computes the sum ``Σ_i M_i V f_i(S)`` for a DEP
    function compute_MM(nep::DEP,S,V)
        Z=-V*S;
        for j=1:size(nep.A,1)
            Z+=nep.A[j]*V*expm(Matrix(-nep.tauv[j]*S))
        end
        return Z
    end

    #  Fetch the Av's, since they are not explicitly stored in DEPs
    function get_Av(nep::DEP)
        local I;
        if (issparse(nep))
            I=speye(eltype(nep.A[1]),nep.n)
        else
            I=eye(eltype(nep.A[1]),nep.n)
        end
        return [I,nep.A...];
    end
    #  Fetch the Fv's, since they are not explicitly stored in DEPs
    function get_fv(nep::DEP)
        fv=Array{Function,1}(size(nep.A,1)+1);
        # First function is -λ
        fv[1] =  S-> -S

        # The other functions are expm(-tauv[j]*S)
        for i=1:size(nep.A,1)
            if (nep.tauv[i]==0) # Zero delay means constant term
                fv[i+1]=  (S-> eye(size(S,1)))
            else
                fv[i+1]=  (S-> expm(-Matrix(nep.tauv[i]*S)))
            end

        end

        return fv;
    end

    ###########################################################
    # Polynomial eigenvalue problem - PEP
    #

    """
    struct PEP <: AbstractSPMF

A polynomial eigenvalue problem (PEP) is defined by the sum the sum ``Σ_i A_i λ^i``, where i = 0,1,2,..., and  all of the matrices are of size n times n.
"""

    struct PEP <: AbstractSPMF
        n::Int
        A::Array   # Monomial coefficients of PEP
    end
"""
    PEP(AA::Array)

Creates a polynomial eigenvalue problem with monomial matrices specified in
AA, which is an array of matrices.

```julia-repl
julia> A0=[1 3; 4 5]; A1=A0+eye(2); A2=ones(2,2);
julia> pep=PEP([A0,A1,A2])
julia> compute_Mder(pep,3)-(A0+A1*3+A2*9)
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0
```
"""
    function PEP(AA::Array)
        n=size(AA[1],1)
        AA=reshape(AA,size(AA,1))
        return PEP(n,AA)
    end

# Computes the sum ``Σ_i M_i V f_i(S)`` for a PEP
    function compute_MM(nep::PEP,S,V)
        if (issparse(nep))
            Z=spzeros(size(V,1),size(V,2))
            Si=speye(size(S,1))
        else
            Z=zeros(size(V))
            Si=eye(size(S,1))
        end
        for i=1:size(nep.A,1)
            Z+=nep.A[i]*V*Si;
            Si=Si*S;
        end
        return Z
    end
    # Use MM to compute Mlincomb for PEPs
    compute_Mlincomb(nep::PEP,λ::Number,
                     V::Union{AbstractMatrix,AbstractVector},a::Vector=ones(size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)

    compute_rf(nep::PEP,x;params...) = compute_rf(Complex128,nep,x;params...)
    function compute_rf{T<:Real}(::Type{T},nep::PEP,x; y=x, target=zero(T), λ0=target,
                           TOL=eps(real(T))*1e3,max_iter=10)

        a=zeros(T,size(nep.A,1))
        for i=1:size(nep.A,1)
            a[i]=dot(y,nep.A[i]*x);
        end
        p=Poly(a);
        rr=real(poly_roots(p));  # Only works if polynomial roots are real
        return rr
    end


# Compute the ith derivative of a PEP
    function compute_Mder(nep::PEP,λ::Number,i::Integer=0)
        if (issparse(nep))
            Z=spzeros(eltype(nep.A[1]),size(nep,1),size(nep,1));
        else
            Z=zeros(eltype(nep.A[1]),size(nep,1),size(nep,1));
        end
        for j=(i+1):size(nep.A,1)
            # Derivatives of monimials
            Z+= nep.A[j]*(λ^(j-i-1)*factorial(j-1)/factorial(j-i-1))
        end
        return Z
    end

    #  Fetch the Av's, since they are not explicitly stored in PEPs
    function get_Av(nep::PEP)
        return nep.A;
    end
    #  Fetch the Fv's, since they are not explicitly stored in PEPs
    function get_fv(nep::PEP)
        fv=Array{Function,1}(size(nep.A,1));
        # Construct monomial functions
        for i=1:size(nep.A,1)
            if (i==1); # optimization for constant and linear term
                fv[1]=(S->eye(size(S,1)));
            elseif (i==2);
                fv[2]=(S->S);
            else
                fv[i]=  S->S^(i-1);
            end
        end
        return fv;
    end


"""
    interpolate([T::DataType=Complex64,] nep::NEP, intpoints::Array)
 Interpolates a NEP in the points intpoints and returns a PEP.\\
 `T` is the DataType in which the PEP should be defined.
"""
    function interpolate(T::DataType, nep::NEP, intpoints::Array)

        n = size(nep, 1)
        d = length(intpoints)

        V = Array{T}(d,d) #Vandermonde matrix
        pwr = ones(d,1)
        for i = 1:d
            V[:,i] = pwr
            pwr = pwr.*intpoints
        end

        if (issparse(nep)) #If Sparse, do elementwise interpolation
            b = Array{SparseMatrixCSC{T},1}(d)
            AA = Array{SparseMatrixCSC{T},1}(d)
            V = factorize(V) # Will be used multiple times, factorize

            for i=1:d
                b[i] = compute_Mder(nep, intpoints[i])
            end

            # OBS: The following lines and hence the  following method assumes that Sparsity-structure is the same!
            nnz_AA = nnz(b[1])
            for i=1:d
                AA[i] = spones(b[1])
            end

            f = zeros(d,1)
            for i = 1:nnz_AA
                for j = 1:d
                    f[j] = b[j].nzval[i]
                end
                a = \(V,f)
                for j = 1:d
                    AA[j].nzval[i] = a[j]
                end
            end

        else # If dense, use Vandermonde
            b = Array{T}(n*d,n)
            AA = Array{Array{T,2}}(d)
            (L, U, p) = lu(V)

            I = speye(n,n)
            LL = kron(L,I)
            UU = kron(U,I)

            for i = 1:d
                b[(1:n)+(i-1)*n,:] =  compute_Mder(nep,intpoints[p[i]])
            end

            A = \(UU, \(LL,b))

            for i = 1:d
                AA[i] = A[(1:n)+(i-1)*n,:]
            end
        end

        return PEP(AA)
    end


    interpolate(nep::NEP, intpoints::Array) = interpolate(Complex128, nep, intpoints)


    """
     interpolate_cheb(nep::NEP,a::Real,b::Real)
  Interpolation in an interval using Chebyshev distribution. Returns a PEP.
  Following Effenberger, Cedric, and Daniel Kressner. "Chebyshev interpolation for nonlinear eigenvalue problems." BIT Numerical Mathematics 52.4 (2012): 933-951.
"""
    function interpolate_cheb(nep::NEP,a::Real,b::Real)
        # Not yet implemented
        # Note: PEP should probably be separated into Mono_PEP and
        # Cheb_PEP depending which inherit from PEP.
        error("Not implemented")
    end




    ###########################################################
    # Rational eigenvalue problem - REP
"""
    struct REP <: AbstractSPMF

A REP represents a rational eigenvalue problem. The REP is defined by the
sum ``Σ_i A_i s_i(λ)/q_i(λ)``, where i = 0,1,2,..., all of the
matrices are of size n times n and s_i and q_i are polynomials.
"""
    struct REP <: AbstractSPMF
        n::Int
        A::Array   # Monomial coefficients of REP
        si::Array  # numerator polynomials
        qi::Array  # demonimator polynomials
    end
    """
    REP(A,poles)

Creates a rational eigenvalue problem. The
constructor takes the matrices A_i and a sequence of poles as input
(not complete).


# Example
```julia-repl
julia> A0=[1 2; 3 4]; A1=[3 4; 5 6];
julia> nep=REP([A0,A1],[1,3]);
julia> compute_Mder(nep,3)
2×2 Array{Float64,2}:
 NaN  NaN
 NaN  NaN
```
"""
    function REP(AA,poles::Array{<:Number,1})

        n=size(AA[1],1)
        AA=reshape(AA,length(AA)) # allow for 1xn matrices
        # numerators
        si=Array{Array{Number,1},1}(length(poles))
        for i =1:size(poles,1)
            si[i]=[1];
        end
        # denominators
        qi=Array{Array{Number,1}}(length(poles))
        for i =1:size(poles,1)
            if poles[i]!=0
                qi[i]=[1,-1/poles[i]];
            else
                qi[i]=[1];
            end
        end
        return REP(n,AA,si,qi)
    end


   function compute_MM(nep::REP,S,V)
        local Z0;
        if (issparse(nep))
            Z=spzeros(size(V,1),size(V,2))
            Si=speye(S)
        else
            Z=zeros(V)
            Si=eye(S)
        end
        # Sum all the elements
        for i=1:size(nep.A,1)
            # compute numerator
            Snum=0*copy(Z);
            Spowj=copy(Si);
            for j=1:size(nep.si[i],1)
                Snum+=Spowj*nep.si[i][j]
                Spowj=Spowj*S;
            end

            # compute denominator
            Sden=0*copy(Z);
            Spowj=copy(Si);
            for j=1:size(nep.qi[i],1)
                Sden+=Spowj*nep.qi[i][j]
                Spowj=Spowj*S;
            end

            # Sum it up
            Z+=nep.A[i]*V*(Sden\Snum)
        end
        return Z
    end
    function compute_Mder(rep::REP,λ::Number,i::Integer=0)
        if (i!=0) # Todo
            error("Higher order derivatives of REP's not implemented")
        end
        S=speye(rep.n)*λ
        V=speye(rep.n);
        return compute_MM(rep,S,V)  # This call can be slow
    end


    #  Fetch the Av's, since they are not explicitly stored in REPs
    function get_Av(nep::REP)
        return nep.A;
    end
    #  Fetch the Fv's, since they are not explicitly stored in REPs
    function get_fv(nep::REP)
        fv=Array{Function,1}(size(nep.qi,1))
        for i=1:size(fv,1)
            fv[i]=S -> (lpolyvalm(nep.qi[i],S)\lpolyvalm(nep.si[i],S))
        end
        return fv
    end

    # Evaluation of matrix polynomial with coefficient a
    function lpolyvalm(a::Array{<:Number,1},S::Array{<:Number,2})
        Sp=eye(S);
        Ssum=zeros(S);
        for j=1:size(a,1)
            Ssum+= a[j]*Sp;
            Sp=Sp*S;
        end
        return Ssum
    end

    function issparse(nep::REP)
        return issparse(nep.A[1])
    end



    #######################################################
    ### Represents a projected NEP
"""
Proj_NEP represents a projected NEP
"""
    abstract type Proj_NEP <: NEP end

"""
    pnep=create_proj_NEP(orgnep::ProjectableNEP)

Create a NEP representing a projected problem. The projection is defined
as the problem ``N(λ)=W^HM(λ)V``
where ``M(λ)`` is represented by `orgnep`. Use
`set_projectionmatrices!()` to specify projection matrices
``V`` and ``W``.
"""
#    function create_proj_NEP(orgnep::ProjectableNEP)
#        if (isa(orgnep,PEP))
#            return Proj_PEP(orgnep);
#        elseif (isa(orgnep,SPMF_NEP))
#            return Proj_SPMF_NEP(orgnep);
#        else
#            error("Projection of this NEP is not available");
#        end
#    end
    function create_proj_NEP(orgnep::ProjectableNEP)
         error("Not implemented. All ProjectableNEP have to implement create_proj_NEP.")
    end
    function create_proj_NEP(orgnep::AbstractSPMF)
         return Proj_SPMF_NEP(orgnep);
    end


    # concrete types for projection of NEPs and PEPs
#    struct Proj_PEP <: Proj_NEP
#        orgnep::PEP
#        V
#        W
#        nep_proj::PEP; # An instance of the projected NEP
#        function Proj_PEP(nep::PEP)
#            this=new(nep);
#            return this
#        end
#    end

    mutable struct Proj_SPMF_NEP <: Proj_NEP
        orgnep::AbstractSPMF
        V
        W
        nep_proj::SPMF_NEP; # An instance of the projected NEP
        orgnep_Av::Array
        orgnep_fv::Array
        function Proj_SPMF_NEP(nep::AbstractSPMF)
            this = new(nep)

            this.orgnep_Av = get_Av(nep)
            if     (size(this.orgnep_Av,1) != 1) && (size(this.orgnep_Av,2) == 1) # Stored as column vector - do nothing
            elseif (size(this.orgnep_Av,1) == 1) && (size(this.orgnep_Av,2) == 1) # It is a single entry - do nothing
            elseif (size(this.orgnep_Av,1) == 1) && (size(this.orgnep_Av,2) != 1) # Stored as a row-vector
                this.orgnep_Av = vec(this.orgnep_Av)
            else
                error("The given array should be a vector but is of size ", size(this.orgnep_Av), ".")
            end

            this.orgnep_fv = get_fv(nep)
            if     (size(this.orgnep_fv,1) != 1) && (size(this.orgnep_fv,2) == 1) # Stored as column vector - do nothing
            elseif (size(this.orgnep_fv,1) == 1) && (size(this.orgnep_fv,2) == 1) # It is a single entry - do nothing
            elseif (size(this.orgnep_fv,1) == 1) && (size(this.orgnep_fv,2) != 1) # Stored as a row-vector
                this.orgnep_fv = vec(this.orgnep_fv)
            else
                error("The given array should be a vector but is of size ", size(this.orgnep_fv), ".")
            end
            return this
        end
    end

"""
    set_projectmatrices!(pnep::Proj_NEP,W,V)
Set the projection matrices for the NEP to W and V, i.e.,
corresponding the NEP: ``N(λ)=W^HM(λ)V``.

# Example:
The following example illustrates that a projection
of a `NEP` is also a `NEP` and we can for instance
call `compute_Mder`on it:
```julia-repl
julia> nep=nep_gallery("pep0")
julia> V=eye(size(nep,1),2);
julia> W=eye(size(nep,1),2);
julia> pnep=create_proj_NEP(nep);
julia> set_projectmatrices!(pnep,W,V);
julia> compute_Mder(pnep,λ)
2×2 Array{Float64,2}:
 -2.03662   13.9777
 -1.35069  -13.0975
julia> compute_Mder(nep,λ)[1:2,1:2]
2×2 Array{Float64,2}:
 -2.03662   13.9777
 -1.35069  -13.0975
```
"""
    function set_projectmatrices!(nep::Proj_SPMF_NEP,W,V)
        ## Sets the left and right projected basis and computes
        ## the underlying projected NEP
        m = size(nep.orgnep_Av,1);
        B = Array{Array{eltype(W),2}}(m);
        for i=1:m
            B[i]=W'*nep.orgnep_Av[i]*V;
        end
        nep.W=W;
        nep.V=V;
        # Keep the sequence of functions for SPMFs
        nep.nep_proj=SPMF_NEP(B,nep.orgnep_fv)
    end


    # Use delagation to the nep_proj
    compute_MM(nep::Union{Proj_SPMF_NEP},par...)=compute_MM(nep.nep_proj,par...)
    # Use MM to compute Mlincomb for SPMFs
    compute_Mlincomb(nep::Proj_SPMF_NEP,λ::Number,
                     V::Union{AbstractMatrix,AbstractVector},a::Vector=ones(size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)
    compute_Mder(nep::Union{Proj_SPMF_NEP},λ::Number)=compute_Mder(nep.nep_proj,λ,0)
    compute_Mder(nep::Union{Proj_SPMF_NEP},λ::Number,i::Integer)=compute_Mder(nep.nep_proj,λ,i)

    function size(nep::Proj_NEP,dim=-1)
        n=size(nep.W,2);
        if (dim==-1)
            return (n,n)
        else
            return n
        end
    end

    function get_Av(nep::Proj_SPMF_NEP)
        return get_Av(nep.nep_proj);
    end

    function get_fv(nep::Proj_SPMF_NEP)
        return orgnep_fv
    end

    function issparse(nep::Proj_NEP)
        return false; # A projected NEP is (essentially) never sparse
    end

    """
        type SumNEP{nep1::NEP,nep2::NEP} <: NEP

SumNEP corresponds to a sum of two NEPs, i.e., if nep is a SumNEP it
is defined by
```math
M(λ)=M_1(λ)+M_2(λ)
```
where M_1 and M_2 are defined by `nep1` and `nep2`.

# Example:
```julia-repl
julia> nep1=DEP([ones(3,3),randn(3,3)])
julia> nep2=PEP([ones(3,3),randn(3,3),randn(3,3)])
julia> sumnep=SumNEP(nep1,nep2);
julia> s=3.0;
julia> M=compute_Mder(sumnep,s);
3×3 Array{Float64,2}:
  8.54014     6.71897   7.12007
 -0.943908  -13.0795   -0.621659
  6.03155    -7.26726  -6.42828
julia> M1=compute_Mder(nep1,s);
julia> M2=compute_Mder(nep2,s);
julia> M1+M2  # Same as M
3×3 Array{Float64,2}:
  8.54014     6.71897   7.12007
 -0.943908  -13.0795   -0.621659
  6.03155    -7.26726  -6.42828
```
"""
    struct GenericSumNEP{NEP1<:NEP,NEP2<:NEP}  <: NEP
        nep1::NEP1
        nep2::NEP2
    end

    struct SPMFSumNEP{NEP1<:AbstractSPMF,NEP2<:AbstractSPMF}  <: AbstractSPMF
        nep1::NEP1
        nep2::NEP2
    end

    AnySumNEP=Union{GenericSumNEP,SPMFSumNEP};
    # Creator functions for SumNEP
    function SumNEP(nep1::AbstractSPMF,nep2::AbstractSPMF)
        return SPMFSumNEP(nep1,nep2);
    end
    function SumNEP(nep1::NEP,nep2::NEP)
        return GenericSumNEP(nep1,nep2);
    end

    # Delegate all interface functions
    size(snep::AnySumNEP)=size(snep.nep1)
    size(snep::AnySumNEP,d)=size(snep.nep1,d)
    compute_Mlincomb(nep::AnySumNEP, λ::Number, V::Union{AbstractMatrix,AbstractVector}) =
        (compute_Mlincomb(nep.nep1, λ, V)+compute_Mlincomb(nep.nep2,λ,V))
    compute_Mder(nep::AnySumNEP, λ::Number,i::Int = 0) =
        (compute_Mder(nep.nep1,λ,i)+compute_Mder(nep.nep2,λ,i))
    compute_MM(nep::AnySumNEP, S::Matrix,V::Matrix) =
        (compute_MM(nep.nep1,S,V)+compute_M(nep.nep2,S,V))

    # For SPMFSumNEP, also delegate the get_Av() and get_fv()
    get_Av(nep::SPMFSumNEP) = [get_Av(nep.nep1); get_Av(nep.nep2)]
    get_fv(nep::SPMFSumNEP) = [get_fv(nep.nep1); get_fv(nep.nep2)]


   #######################################################
   ### Functions in common for many NEPs in NEPTypes

   #
"""
    size(nep::NEP,dim=-1)
 Overloads the size functions for NEPs storing size in nep.n
"""
    function size(nep::Union{DEP,PEP,REP,SPMF_NEP},dim)
        return nep.n
    end


    function size(nep::Union{DEP,PEP,REP,SPMF_NEP})
        return (nep.n,nep.n)
    end

"""
    issparse(nep)
Returns true/false if the NEP is sparse (if compute_Mder() returns sparse)
"""
    function issparse(nep::Union{DEP,PEP,REP,SPMF_NEP})
        return issparse(nep.A[1])
    end


    function get_Av(nep::Union{SPMF_NEP,PEP,REP})
        return nep.A;
    end
    function get_fv(nep::SPMF_NEP)
        return nep.fi;
    end


    include("nep_transformations.jl")

    # structure exploitation for DEP (TODO: document this)
    function compute_Mlincomb(nep::DEP,λ::T,V::Matrix{T},
                              a::Vector{T}=ones(T,size(V,2))) where {T<:Number}
        n=size(V,1); k=1
        try k=size(V,2) end
        Av=get_Av(nep)
        broadcast!(*,V,V,a.')
        z=zeros(T,n)
        for j=1:length(nep.tauv)
            w=Array{T,1}(exp(-λ*nep.tauv[j])*(-nep.tauv[j]).^(0:k-1))
            z[:]+=Av[j+1]*(V*w);
        end
        if k>1 z[:]-=view(V,:,2:2) end
        z[:]-=λ*view(V,:,1:1);
        return z
    end
    # Automatically promote to complex if λ is real
    function compute_Mlincomb(nep::DEP,λ::T,V::Array{Complex{T},2},a::Vector{Complex{T}}=ones(Complex{T},size(V,2))) where T<:Real
        return compute_Mlincomb(nep,complex(λ),V,a)
    end
    # Allow vector-valued V
    function compute_Mlincomb(nep::DEP,λ::Number,V::Vector{T},a::Vector{T}=ones(T,1)) where T<:Number
        return compute_Mlincomb(nep,λ,reshape(V,size(V,1),1),a)
    end
end
