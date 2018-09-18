module NEPTypes
    using ..NEPCore
    using SparseArrays
    using LinearAlgebra
    using PolynomialZeros
    using Polynomials
    using InteractiveUtils

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
    export get_Av
    export get_fv

    export interpolate
    export set_projectmatrices!



    # We overload these
    import ..NEPCore.compute_Mder
    import ..NEPCore.compute_Mlincomb
    import ..NEPCore.compute_Mlincomb!
    import ..NEPCore.compute_MM
    import ..NEPCore.compute_resnorm
    import ..NEPCore.compute_rf
    import Base.size
    import SparseArrays.issparse



    #
    """
    abstract ProjectableNEP <: NEP

A ProjectableNEP is a NEP which can be projected, i.e., one can construct the problem W'*M(λ)Vw=0 with the Proj_NEP. A NEP which is of this must have the function `create_proj_NEP(orgnep::ProjectableNEP)` implemented. This function must return a `Proj_NEP` See also `set_projectmatrices!".

# Example:
julia> nep=nep_gallery("dep0");
julia> typeof(nep)
DEP{Float64,Array{Float64,2}}
julia> isa(nep,ProjectableNEP)
true
julia> projnep=create_proj_NEP(nep);
julia> e1 = Matrix(1.0*I,size(nep,1),1);
julia> set_projectmatrices!(projnep,e1,e1);
julia> compute_Mder(nep,3.0)[1,1]
-2.315345215259418
julia> compute_Mder(projnep,3.0)
1×1 Array{Float64,2}:
 -2.315345215259418
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
         As::Vector{SparseMatrixCSC{<:Number,Int}}  # 'A' matrices with sparsity pattern of all matrices combined

         check_consistency::Bool # input checking
    end

    SPMF_NEP(n, A, fi, Schur_factorize_before, Zero,set_check_consistency) =
        SPMF_NEP(n, A, fi, Schur_factorize_before, Zero, Vector{SparseMatrixCSC{Float64,Int}}(),check_consistency)

"""
     SPMF_NEP(AA, fii, Schur_fact = false, use_sparsity_pattern = true, check_consistency=true)

Creates a `SPMF_NEP` consisting of matrices `AA` and functions `fii`. `fii` must be an array of functions defined for matrices and numbers. `AA` is an array of matrices. `Schur_fact` specifies if the computation of `compute_MM` should be done by first pre-computing a Schur-factorization (which can be faster). If `use_sparsity_pattern` is true, and the `AA` matrices are sparse, each matrix will be stored with a sparsity pattern matching the union of all `AA` matrices. This leads to more efficient calculation of `compute_Mder`. If the sparsity patterns are completely or mostly distinct, it may be more efficient to set this flag to false. If check_consistency is true the input
checking will be performed.

```julia-repl
julia> A0=[1 3; 4 5]; A1=[3 4; 5 6];
julia> id_op=S -> one(S)
julia> exp_op=S -> exp(S)
julia> nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
julia> compute_Mder(nep,1)-(A0+A1*exp(1))
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0
```
"""
     function SPMF_NEP(AA::Vector{<:AbstractMatrix}, fii::Vector{<:Function},
            Schur_fact = false, use_sparsity_pattern = true, check_consistency=true)

            T=Float64;
            if (check_consistency)
                for t=1:length(fii)
                    # Scalar leads to scalars:
                    s=one(T);
                    ci=@code_typed(fii[t](s)) # ci[end] gives the return type
                    if (ci[end] <: Number)
                        # nothing
                    else
                        @warn "It seems you have provided valid matrix-functions for defining SPMF_NEP. The functions fii should return a scalar if evaluated in a scalar and a matrix if evaluated in a matrix. If you want to disable to input checking, set check_consistency=false in SPMF_NEP."
                        #error("The given function does not return a scalar if evaluated in a scalar")
                    end
                    S=ones(T,2,2);
                    ci=@code_typed(fii[t](S))
                    if (ci[end] <: Matrix)
                        # nothing
                    else
                        @warn "It seems you have provided valid matrix-functions for
                        defining SPMF_NEP. The functions fii should return a scalar if evaluated in a scalar and a matrix if evaluated in a matrix. If you want to disable to input checking, set check_consistency=false in SPMF_NEP."
                        #error("The given function does not return a matrix if evaluated in a matrix")
                    end
                end
            end

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


         T = eltype(AA[1])
         As = Vector{SparseMatrixCSC{T,Int}}()

         if use_sparsity_pattern && issparse(AA[1])
             # Merge the sparsity pattern of all matrices without dropping any zeros
             Zero = LinearAlgebra.fillstored!(copy(AA[1]), 1)
             for i = 2:size(AA,1)
                 Zero += LinearAlgebra.fillstored!(copy(AA[i]), 1)
             end
             Zero = T.(Zero)
             Zero.nzval[:] .= T(0)

             # Create a copy of each matrix with the sparsity pattern of all matrices combined
             @inbounds for A in AA
                 S = copy(Zero)

                 for col = 1:size(A, 2)
                    for j in nzrange(A, col)
                        S[A.rowval[j], col] = A.nzval[j]
                    end
                 end

                 push!(As, S)
             end
         else
             Zero=zeros(n,n)
         end


         this=SPMF_NEP(n,AA,fii,Schur_fact,Zero,As,true);
         return this
    end
    function SPMF_NEP(n) # Create an empty NEP of size n x n
         Z=zeros(n,n)
         return SPMF_NEP(n,Vector{Matrix}(),Vector{Function}(),false,Z,false);
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
                if (norm(Sd .- Sd[1])==0) # Optimize further if S is a multiple of identity
                    Fid=nep.fi[i](reshape([Sd[1]],1,1))[1]*ones(size(Sd,1))
                else  # Diagonal but not constant
                    Fid=zeros(ComplexF64,size(S,1))
                    for j=1:size(S,1)
                        Fid[j]=nep.fi[i](reshape([Sd[j]],1,1))[1]
                    end
                end
                Fi=sparse(Diagonal(Fid))
            else  # Otherwise just compute the matrix function operation
                if(nep.Schur_factorize_before)
                    Fi= Q*nep.fi[i](T)*Q'
                else
                    Fi=nep.fi[i](S)
                end
            end
            ## Sum it up


            VFi=V*Fi;
            if (isa(nep.A[i],SubArray) && (eltype(nep.A[i]) != eltype(VFi)))
                # 2018-09-14: SubArray x Matrix of different types do not work? Check in new version of Julia/base at later point (TODO)
                # https://discourse.julialang.org/t/subarray-x-matrix-multipliciation-of-different-eltype-fails/14950
                # Workaround by making a matrix copy.
                Z=Z+copy(nep.A[i])*VFi;
            else
                Z=Z+nep.A[i]*VFi;
            end

        end
        return Z
    end


    function compute_Mder(nep::SPMF_NEP,λ::Number)

        # figure out the return type, as the greatest type of all input
        x = map(i -> nep.fi[i](reshape([λ],1,1))[1], 1:length(nep.fi))
        Tx = mapreduce(eltype, promote_type, x)
        TA=mapreduce(eltype, promote_type, nep.A); # Greatest type of all A-matrices
        TZ=promote_type(TA,Tx)  # output type
        if isempty(nep.As)
            # Full matrices
            Z = zeros(TZ,size(nep,1),size(nep,1));
            for k=1:size(nep.A,1)
                Z += nep.A[k] * x[k]
            end
            return Z
        else
            Z = SparseMatrixCSC(nep.As[1].m, nep.As[1].n,
                                nep.As[1].colptr, nep.As[1].rowval,
                                convert.(TZ, nep.As[1].nzval .* x[1]))
            for k = 2:length(nep.As)
                Z.nzval .+= nep.As[k].nzval .* x[k]
            end

            return Z
        end
    end

    function compute_Mder(nep::SPMF_NEP,λ::Number,i::Integer)
        if (i==0)
            return compute_Mder(nep,λ);
        else
            # This is typically slow for i>1 (can be optimized by
            # treating the SPMF-terms individually)
            return compute_Mder_from_MM(nep,λ,i)
        end
    end


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

# Example:
julia> A0=randn(3,3); A1=randn(3,3);
julia> tauv=[0,0.2] # Vector with delays
julia> dep=DEP([A0,A1],tauv)
julia> λ=3.0;
julia> M1=compute_Mder(dep,λ)
julia> M2=-λ*I+A0+A1*exp(-tauv[2]*λ)
julia> norm(M1-M2)
0.0
"""
    struct DEP{Z<:Real, T<:AbstractMatrix} <: AbstractSPMF
        n::Int
        A::Array{T,1}     # An array of matrices (full or sparse matrices)
        tauv::Vector{Z}   # the delays (which are always real)
    end
    function DEP(AA::Vector{T},tauv::Vector=[0,1.0]) where {T<:AbstractMatrix}
        n=size(AA[1],1)
        if (!isreal(tauv))
            error("Incorrect construction of DEP. The delays need to be real.")
        end

        # Note: this enforces that eltype(tauv)=real(eltype(A[1]))
        tauvconv::Vector{real(eltype(AA[1]))}=Vector{real(eltype(AA[1]))}(tauv);

        this=DEP{eltype(tauvconv),T}(n,AA,tauvconv);
        return this;
    end

# Compute the ith derivative of a DEP
    function compute_Mder(nep::DEP,λ::Number,i::Integer=0)
        local M, J
        # T is eltype(nep.A[1]) unless λ complex, then T is complex(eltype(nep.A[1]))
        # T can be determined compile time, since DEP parametric type
        T=isa(λ,Complex) ? complex(eltype(nep.A[1])) : eltype(nep.A[1]);

        TM=promote_type(eltype(λ),T); # output type
        if (issparse(nep.A[1])) # Can be determined compiled time since DEP parametric type
            M=spzeros(TM,nep.n,nep.n)
            J=sparse(TM(1)I, nep.n, nep.n)
        else
            M=zeros(TM,nep.n,nep.n)
            J=Matrix{TM}(I, nep.n, nep.n)
        end
        if i==0; M=-λ*J;  end
        if i==1; M=-one(TM)*J; end
        for j=1:size(nep.A,1)
            a=exp(-nep.tauv[j]*λ)*(-nep.tauv[j])^i;
            M += nep.A[j]*a
        end
        return M
    end




#  Computes the sum ``Σ_i M_i V f_i(S)`` for a DEP
    function compute_MM(nep::DEP,S,V)

        T1=promote_type(promote_type(eltype(S),eltype(V)),eltype(nep.A[1]));
        T=promote_type(T1,eltype(nep.tauv));

        Z::Matrix{T}=-V*S;
        for j=1:size(nep.A,1)
            Z+=nep.A[j]*V*exp(Matrix(-nep.tauv[j]*S))
        end
        return Z
    end

    #  Fetch the Av's, since they are not explicitly stored in DEPs
    function get_Av(nep::DEP)
        local J
        if issparse(nep)
            J = sparse(eltype(nep.A[1])(1)I, nep.n, nep.n)
        else
            J = Matrix{eltype(nep.A[1])}(I, nep.n, nep.n)
        end
        return [J, nep.A...]
    end
    #  Fetch the Fv's, since they are not explicitly stored in DEPs
    function get_fv(nep::DEP)
        fv = Array{Function,1}(undef, size(nep.A,1)+1)
        # First function is -λ
        fv[1] = S -> -S

        # The other functions are exp(-tauv[j]*S)
        for i=1:size(nep.A,1)
            if nep.tauv[i] == 0 # Zero delay means constant term
                fv[i+1] = S -> Matrix(1.0I, size(S,1), size(S,1))
            else
                fv[i+1] = S -> exp(-Matrix(nep.tauv[i] * S))
            end
        end

        return fv
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
julia> A0=[1 3; 4 5]; A1=A0.+one(2); A2=ones(2,2);
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
        T=promote_type(promote_type(eltype(nep.A[1]),eltype(S)),eltype(V))
        local Z
        if (issparse(nep))
            Z=spzeros(T,size(V,1),size(V,2))
            Si=sparse(one(T)*I, size(S,1), size(S,1))
        else
            Z=zeros(T,size(V,1),size(V,2))
            Si=Matrix(one(T)*I, size(S,1), size(S,1))
        end
        for i=1:size(nep.A,1)
            Z+=nep.A[i]*V*Si;
            Si=Si*S;
        end
        return Z
    end
    # Use MM to compute Mlincomb for PEPs
    compute_Mlincomb(nep::PEP,λ::Number,
                     V::Union{AbstractMatrix,AbstractVector},a::Vector=ones(eltype(V),size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)

    compute_rf(nep::PEP,x;params...) = compute_rf(ComplexF64,nep,x;params...)
    function compute_rf(::Type{T},nep::PEP,x; y=x, target=zero(T), λ0=target,
                        TOL=eps(real(T))*1e3,max_iter=10) where T<:Real

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
        fv=Vector{Function}(undef, size(nep.A,1))
        # Construct monomial functions
        for i=1:size(nep.A,1)
            if (i==1); # optimization for constant and linear term
                fv[1] = S -> Matrix(1.0I, size(S, 1), size(S, 1))
            elseif (i==2);
                fv[2]=(S->S);
            else
                fv[i]=  S->S^(i-1);
            end
        end
        return fv;
    end


"""
    interpolate([T=ComplexF64,] nep::NEP, intpoints::Array)
 Interpolates a NEP in the points `intpoints` and returns a `PEP`.\\
 `T` is the DataType in which the matrices of the PEP should be defined.
"""
    interpolate(nep::NEP, intpoints::Array) = interpolate(ComplexF64, nep, intpoints)
    function interpolate(::Type{T}, nep::NEP, intpoints::Array) where {T<:Number}

        n = size(nep, 1)
        d = length(intpoints)

        V = zeros(T,d,d) #Vandermonde matrix
        pwr = ones(d,1)
        for i = 1:d
            V[:,i] = pwr
            pwr = pwr.*intpoints
        end

        if (issparse(nep)) #If Sparse, do elementwise interpolation
            b = Vector{SparseMatrixCSC{T}}(undef, d)
            AA = Vector{SparseMatrixCSC{T}}(undef, d)
            V = factorize(V) # Will be used multiple times, factorize

            for i=1:d
                b[i] = compute_Mder(nep, intpoints[i])
            end

            # OBS: The following lines and hence the following method assumes that Sparsity-structure is the same!
            nnz_AA = nnz(b[1])
            for i=1:d
                AA[i] = copy(b[1])
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
            b = zeros(T,n*d,n)
            AA = Vector{Matrix{T}}(undef,d)
            (L, U, p) = lu(V)

            LL = kron(L, SparseMatrixCSC(I,(n,n)))
            UU = kron(U, SparseMatrixCSC(I,(n,n)))

            for i = 1:d
                b[(1:n).+(i-1)*n,:] =  compute_Mder(nep,intpoints[p[i]])
            end

            A = \(UU, \(LL,b))

            for i = 1:d
                AA[i] = A[(1:n).+(i-1)*n,:]
            end
        end

        return PEP(AA)
    end
  # TODO: Implement interpolation similar to Effenberger and Kressner. "Chebyshev interpolation for nonlinear eigenvalue problems." BIT Numerical Mathematics 52.4 (2012): 933-951.




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
        si = Vector{Vector{Number}}(undef, length(poles))
        for i =1:size(poles,1)
            si[i]=[1];
        end
        # denominators
        qi = Vector{Vector{Number}}(undef, length(poles))
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
            Si = SparseMatrixCSC{eltype(S)}(I, size(S))
        else
            Z = zero(V)
            Si = Matrix{eltype(S)}(I, size(S))
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
        S = SparseMatrixCSC(λ*I, rep.n, rep.n)
        V = SparseMatrixCSC(1.0I, rep.n, rep.n)
        return compute_MM(rep,S,V)  # This call can be slow
    end


    #  Fetch the Av's, since they are not explicitly stored in REPs
    function get_Av(nep::REP)
        return nep.A;
    end
    #  Fetch the Fv's, since they are not explicitly stored in REPs
    function get_fv(nep::REP)
        fv = Vector{Function}(undef, size(nep.qi, 1))
        for i=1:size(fv,1)
            fv[i]=S -> (lpolyvalm(nep.qi[i],S)\lpolyvalm(nep.si[i],S))
        end
        return fv
    end

    # Evaluation of matrix polynomial with coefficient a
    function lpolyvalm(a::Array{<:Number,1},S::Array{<:Number,2})
        Sp = Matrix{eltype(S)}(I, size(S))
        Ssum = zero(S)
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
    pnep=create_proj_NEP(orgnep::ProjectableNEP[,maxsize [,T]])

Create a NEP representing a projected problem. The projection is defined
as the problem ``N(λ)=W^HM(λ)V`` where ``M(λ)`` is represented by `orgnep`.
The optional parameter `maxsize` determines how large the projected
problem can be and `T` determines which Number type to use (default `ComplexF64`).
These are needed for memory allocation reasons.
Use `set_projectmatrices!()` to specify projection matrices
``V`` and ``W``.
"""
    function create_proj_NEP(orgnep::ProjectableNEP)
         error("Not implemented. All ProjectableNEP have to implement create_proj_NEP.")
    end
    #    function create_proj_NEP(orgnep::ProjectableNEP)
    #        if (isa(orgnep,PEP))
    #            return Proj_PEP(orgnep);
    #        elseif (isa(orgnep,SPMF_NEP))
    #            return Proj_SPMF_NEP(orgnep);
    #        else
    #            error("Projection of this NEP is not available");
    #        end
    #    end
    function create_proj_NEP(orgnep::AbstractSPMF,
                             maxsize::Int=min(size(orgnep,1),
                                              max(round(Int,size(orgnep,1)/10),10)),
                             T::Type{<:Number}=ComplexF64)
        return Proj_SPMF_NEP(orgnep,maxsize,T);
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
        nep_proj::SPMF_NEP; # An instance of the projected NEP
        orgnep_Av::Vector
        orgnep_fv::Vector
        projnep_B_mem::Vector # A vector of matrices
        function Proj_SPMF_NEP(nep::AbstractSPMF,maxsize::Int,T)
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

            this.projnep_B_mem=Vector{Matrix{T}}(undef,size(this.orgnep_fv,1));
            for k=1:size(this.orgnep_fv,1)
                this.projnep_B_mem[k]=zeros(T,maxsize,maxsize);
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
julia> V=Matrix(1.0*I,size(nep,1),2);
julia> W=Matrix(1.0*I,size(nep,1),2);
julia> pnep=create_proj_NEP(nep);
julia> set_projectmatrices!(pnep,W,V);
julia> compute_Mder(pnep,3.0)
2×2 Array{Float64,2}:
 -2.03662   13.9777
 -1.35069  -13.0975
julia> compute_Mder(nep,3.0)[1:2,1:2]
2×2 Array{Float64,2}:
 -2.03662   13.9777
 -1.35069  -13.0975
```
"""
    function set_projectmatrices!(nep::Proj_SPMF_NEP,W,V)
        ## Sets the left and right projected basis and computes
        ## the underlying projected NEP
        m = size(nep.orgnep_Av,1);
        T=eltype(eltype(nep.projnep_B_mem));
        k=size(V,2);
        # Compute first matrix beforhand to determine type
        B1=view(Matrix{T}(copy(W')*nep.orgnep_Av[1]*V),1:k,1:k);
        T_sub = typeof(B1)
        # The coeff matrices for the SPMF_NEP created in the end
        B = Vector{T_sub}(undef,m);
        B[1]=B1;
        for i=2:m # From 2 since we already computed the first above
            nep.projnep_B_mem[i][1:k,1:k]=copy(W')*nep.orgnep_Av[i]*V;
            B[i]=view(nep.projnep_B_mem[i],1:k,1:k);
        end
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

    function size(nep::Proj_NEP,dim)
        n = size(nep.nep_proj,2);
        return n
    end
    function size(nep::Proj_NEP)
        n = size(nep.nep_proj,2);
        return (n,n)
    end

    function get_Av(nep::Proj_SPMF_NEP)
        return get_Av(nep.nep_proj);
    end

    function get_fv(nep::Proj_SPMF_NEP)
        return nep.orgnep_fv
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
        (compute_MM(nep.nep1,S,V)+compute_MM(nep.nep2,S,V))

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
    function compute_Mlincomb(nep::DEP,λ::Number,V::Union{AbstractMatrix,AbstractVector},
                              a::Vector=ones(eltype(V),size(V,2)))
        n=size(V,1); k=size(V,2);
        Av=get_Av(nep)
        # determine type froom greates type of (eltype probtype and λ)
        T=promote_type(promote_type(eltype(V),typeof(λ)),eltype(Av[1]))
        if (k>1)
            broadcast!(*,V,V,transpose(a))
        else
            V=V*a[1];
        end


        z=zeros(T,n)
        for j=1:length(nep.tauv)
            w=Array{T,1}(exp(-λ*nep.tauv[j])*(-nep.tauv[j]) .^(0:k-1))
            if k>1
                z[:]+=Av[j+1]*(V*w);
            else
                z[:]+=Av[j+1]*(V*w[1]);
            end
        end
        if k>1 z[:]-=view(V,:,2:2) end
        z[:]-=λ*view(V,:,1:1);
        return z
    end

    function compute_Mlincomb(
                        nep::SPMF_NEP,
                        λ::Number,
                        V::Union{AbstractMatrix,AbstractVector},
                        a::Vector=ones(size(V,2)))


        if ndims(V)==2
            n,k=size(V);
        else
            n=size(V,1); k=1;
        end

    	# we need to assume that the elements of a are different than zero.
    	V[:,findall(x->x==0,a)] .= 0
    	a[findall(x->x==0,a)] .= 1
    	S=diagm(0 => λ*ones(eltype(V),k)) + diagm(-1 => (a[2:k]./a[1:k-1]).*(1:k-1))

        z=zeros(eltype(V),n)
        if ndims(V)==1
            for i=1:size(nep.A,1)
                Fi=nep.fi[i](S);
                z=z .+ nep.A[i]*(V*Fi);
            end
        else
            for i=1:size(nep.A,1)
                Fi=nep.fi[i](S)[:,1];
                z=z .+ nep.A[i]*(V*Fi);
            end
    	end

    	return a[1]*reshape(z,size(z,1))
    end



end
