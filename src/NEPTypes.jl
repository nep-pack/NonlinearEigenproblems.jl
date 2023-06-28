module NEPTypes
    using ..NEPCore
    using SparseArrays
    using LinearAlgebra
    using InteractiveUtils

    # Specializalized NEPs
    export ProjectableNEP
    export DEP
    export PEP
    export REP
    export SPMF_NEP

    export AbstractSPMF
    export SumNEP, SPMFSumNEP, GenericSumNEP
    export Proj_NEP
    export Proj_SPMF_NEP
    export DerSPMF

    export create_proj_NEP
    export get_Av
    export get_fv

    export interpolate
    export set_projectmatrices!
    export expand_projectmatrices!

    export compute_rf


    # We overload these
    import ..NEPCore.compute_Mder
    import ..NEPCore.compute_Mlincomb
    import ..NEPCore.compute_Mlincomb!
    import ..NEPCore.compute_MM
    import ..NEPCore.compute_resnorm
    import Base.size
    import SparseArrays.issparse


    # This is workaround for Julia bug https://github.com/JuliaLang/julia/issues/29224
    # This can be removed if the fix is backported in julia 1.0.1
    import Base.*
    function (*)(A::SubArray{<:Complex},B::Matrix{<:Real})
        return copy(A)*B;
    end

    "Returns the greatest type of all array elements."
    promote_typeof(A::AbstractArray) = isempty(A) ? eltype(A) : mapreduce(typeof, promote_type, A)

    "Returns the greatest element type of all array elements."
    promote_eltype(A::AbstractArray{<:AbstractArray}) = isempty(A) ? eltype(eltype(A)) : mapreduce(eltype, promote_type, A)

    #
"""
    abstract ProjectableNEP <: NEP

A ProjectableNEP is a NEP which can be projected, i.e., one can construct the problem ``W'*M(λ)Vw=0`` with the [`Proj_NEP`](@ref).
A NEP which is of this must have the function [`create_proj_NEP(orgnep::ProjectableNEP)`](@ref) implemented.
This function must return a `Proj_NEP`.

 See also: [`set_projectmatrices!`](@ref).

# Example:
```julia-repl
julia> nep=nep_gallery("dep0");
julia> typeof(nep)
DEP{Float64,Array{Float64,2}}
julia> isa(nep,ProjectableNEP)
true
julia> projnep=create_proj_NEP(nep);
julia> e1 = Matrix(1.0*I,size(nep,1),1);
julia> set_projectmatrices!(projnep,e1,e1);
julia> compute_Mder(nep,3.0)[1,1]
-2.942777908030041
julia> compute_Mder(projnep,3.0)
1×1 Array{Complex{Float64},2}:
 -2.942777908030041 + 0.0im
```
"""
    abstract type ProjectableNEP <: NEP end

    #######################################################
    ### Sum of products matrices and functions

"""
    abstract  AbstractSPMF <: ProjectableNEP

An AbstractSPMF is an abstract class representing NEPs which can be represented
as a sum of products of matrices and functions ``M(λ)=Σ_i A_i f_i(λ)``,
where i = 0,1,2,..., all of the matrices are of size ``n×n`` and ``f_i`` are functions.

Any AbstractSPMF has to have implementations of [`get_Av()`](@ref) and [`get_fv()`](@ref) which return the
functions and matrices.
"""
    abstract  type AbstractSPMF{T} <: ProjectableNEP end # See issue #17

"""
    get_Av(nep::AbstractSPMF)

Returns an array of matrices ``A_i`` in the AbstractSPMF: ``M(λ)=Σ_i A_i f_i(λ)``
"""
    function get_Av(nep::AbstractSPMF) # Dummy function which enforces that you have to implement
        error("You need to implement get_Av for all AbstractSPMFs")
    end
"""
    get_Av(nep::AbstractSPMF)

Returns an Array of functions (that can be evaluated both as scalar and matrix functions) ``f_i`` in the AbstractSPMF: ``M(λ)=Σ_i A_i f_i(λ)``
"""
    function get_fv(nep::AbstractSPMF) # Dummy function which enforces that you have to implement
        error("You need to implement get_fv for all AbstractSPMFs")
    end

"""
    struct SPMF_NEP{T<:AbstractMatrix,Ftype}  <: AbstractSPMF{T}
    function SPMF_NEP(AA, fii [,check_consistency=true] [,Schur_fact = false]
                      [,align_sparsity_patterns = false] [,Ftype=ComplexF64])

An `SPMF_NEP` is a `NEP` defined by a *S*um of *P*roducts of
*M*atrices and *F*unctions,
i.e.,
```math
M(λ)=∑_i A_i f_i(λ).
```
All of the matrices ``A_0,...`` are of size ``n×n``
and ``f_i`` are a functions. The  functions ``f_i`` must be defined
for matrices in the standard matrix function sense.
The constructor creates a `SPMF_NEP` consisting
of matrices `AA` and functions `fii`.

# Parameters

* `AA` is a `Vector` of matrices. The matrices have to be of the same type. If you need a NEP with different types you can use [`SumNEP`](@ref) to construct a sum of two `SPMF_NEP`.

* `fii` is a `Vector` of functions. Each function takes one parameter `S`. The functions must be available both as a scalar valid function and a matrix function. If `S` is a square matrix, `fii[k](S)` musst also be a square matrix. If `S` is a scalar `fii[k](S)` is a scalar.

* `check_consistency` (default `true`) determines if we should initiate by running tests to verify that the `fii` satisfies the conditions that every function is valid both for matrices and scalars. This is done by using `@code_typed` and the functions need to be type-stable in that sense.

* `align_sparsity_patterns` (default `false`) has effect only for sparse matrices (`SparseMatrixCSC`). If `align_sparsity_patterns=true` the `SparseMatrixCSC` matrices will be replaced by equivalent `SparseMatrixCSC` matrices where the `colptr` and `rowval` are identical. This increases the speed of some functions, e.g., `compute_Mder`. If `align_sparsity_patterns=true` the matrices in the NEP should be considered read only. If the sparsity patterns are completely or mostly distinct, it may be more efficient to set this flag to false.

* `Ftype` (default `ComplexF64`) determines an underlying type of the functions. The output of any function should be "smaller" than the promoted type of the input and `Ftype`. More precisely, if `F=fii[k]`, then the type logic is as follows `eltype(F(λ))=promote_type(eltype(λ),Ftype)`.

* `Schur_fact` (default `false`) determines if the `compute_MM` function should triangularize the matrix before carrying out the computation. This can be faster for large matrices.




# Example
```julia-repl
julia> A0=[1 3; 4 5]; A1=[3 4; 5 6];
julia> id_op=S -> one(S) # Note: We use one(S) to be valid both for matrices and scalars
julia> exp_op=S -> exp(S)
julia> nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
julia> compute_Mder(nep,1)-(A0+A1*exp(1))
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0
```

"""
    struct SPMF_NEP{T<:AbstractMatrix,Ftype}  <: AbstractSPMF{T}
    # Logic behind Ftype:
    #  The eltype(F(λ))=promote_type(eltype(λ),Ftype)
        n::Int
        A::Vector{T}   # Array of Array of matrices
        fi::Vector{Function}  # Array of functions
        Schur_factorize_before::Bool # Tells if you want to do the Schur-factorization
        sparsity_patterns_aligned:: Bool
    end



    # Alias for SPMFs which contain sparse matrices
    SPMF_NEP_sparse{T<:AbstractSparseMatrix,Ftype}=SPMF_NEP{T,Ftype}


    function SPMF_NEP(AA::Vector{<:AbstractMatrix}, fii::Vector{<:Function};
                       Schur_fact = false,
                       check_consistency=true, Ftype=ComplexF64,
                       align_sparsity_patterns=false)

         matrices_are_sparse = issparse.(AA)
         if !(all(matrices_are_sparse) || all(.!matrices_are_sparse))
             error("Mixing sparse and dense matrices is not allowed in SPMF_NEP. Either use a consistent " *
                "format, or split your problem into a sparse and a dense part and use SumNEP.")
         end

         T=Float64;
         if (check_consistency)
             for t=1:length(fii)
                 # Scalar leads to scalars:
                 s=one(T);
                 ci=@code_typed(fii[t](s)) # ci[end] gives the return type
                 if !(ci[end] <: Number)
                     @warn "It seems you have not provided valid matrix-functions for defining SPMF_NEP. The functions fii should return a scalar if evaluated in a scalar and a matrix if evaluated in a matrix. If you want to disable to input checking, set check_consistency=false in SPMF_NEP."
                 end
                 S=ones(T,2,2);
                 ci=@code_typed(fii[t](S))
                 if !(ci[end] <: Matrix)
                     @warn "It seems you have not provided valid matrix-functions for defining SPMF_NEP. The functions fii should return a scalar if evaluated in a scalar and a matrix if evaluated in a matrix. If you want to disable to input checking, set check_consistency=false in SPMF_NEP."
                 end
             end
         end

         if (size(AA,1)==0)
             return SPMF_NEP(0); # Create empty SPMF_NEP.
         end
         n=size(AA[1],1);


         if(size(AA,1) != size(fii,1))
             error("Inconsistency: Number of supplied matrices = ", length(AA), " but the number of supplied functions are = ", length(fii))
         end
         for i = 2:size(AA,1)
             if size(AA[i]) != size(AA[1])
                 error("The dimensions of the matrices mismatch: size(AA[1])",
                       size(AA[1]),"!=",size(AA[i]),"=size(AA[",i,"])")
             end
         end

         if !(eltype(AA) <: SparseMatrixCSC)
             # Dense
             this=SPMF_NEP{promote_typeof(AA),Ftype}(n,AA,fii,Schur_fact,false);
         else
             # Sparse: Potentially do the joint sparsity pattern trick.
             if (!align_sparsity_patterns)
                 # No aligning, just create it
                 this=SPMF_NEP{promote_typeof(AA),Ftype}(n,AA,fii,Schur_fact,false);
             else
                 TT = eltype(AA[1])
                 As=form_aligned_sparsity_patterns(AA,TT)
                 this=SPMF_NEP{promote_typeof(As),Ftype}(n,As,fii,Schur_fact,true);
             end
         end
         return this
    end
    function SPMF_NEP(n) # Create an empty NEP of size n x n
        Z=zeros(n,n)
        return SPMF_NEP{AbstractMatrix,Complex}(n,Vector{Matrix}(),Vector{Function}(),false);
    end

""" Return a vector of sparse matrices which are the same as AA but also have the same sparsity pattern. """
    function form_aligned_sparsity_patterns(AA::Vector{<:SparseMatrixCSC},T)

        # Create a Zero matrix with the correct sparsity pattern.
        Zero = LinearAlgebra.fillstored!(copy(AA[1]), 1)
        for i = 2:size(AA,1)
            Zero .+= LinearAlgebra.fillstored!(copy(AA[i]), 1)
        end
        Zero = T.(Zero)
        Zero.nzval[:] .= T(0) # Set the nonzero elements to zero (keep pattern)


        # Create an vector of sparse matrices. Set the elements corresponding to
        # AA[i] in each of them, such that As[i]=AA[i] but with different sparsity patterns.
        As = Vector{SparseMatrixCSC{T,Int}}()
        @inbounds for A in AA
            # Create a new zero matrix. Note all S-matrices have the same colptr, so
            # modification of As[i] will change all. The As[i] matrices should therefore
            # be viewed as read-only objects.
            S = SparseMatrixCSC(Zero.m, Zero.n,
                                Zero.colptr, Zero.rowval,
                                fill(zero(T),size(Zero.nzval,1)))
            # Set the nonzero elements of A in S
            for col = 1:size(A, 2)
                for j in nzrange(A, col)
                    S[A.rowval[j], col] = A.nzval[j]
                end
            end
            push!(As, S)
        end
        return As;
    end

    function compute_MM(nep::SPMF_NEP{T,Ftype},S,V) where {T,Ftype}

        AA=get_Av(nep);
        ff=get_fv(nep);
        m=size(ff,1); n=size(nep,1); p=size(S,1);

        # Type logic including Ftype
        FStype=promote_type(eltype(S),Ftype) # eltype of f(S)
        T0=promote_type(eltype(V),eltype(AA[1]),FStype) # Output type


        # Always return a dense matrix
        Z=zeros(T0,n,p);

        # Precompute schur factorization (temporarily disabled)
        #if (nep.Schur_factorize_before) #Optimize if allowed to factorize before
        #    (TT, Q, ) = schur(S)  # Currently not used
        #end

        VFi=Matrix{promote_type(FStype,eltype(V))}(undef,n,p); # Type of V*f(S)
        local Fi=Matrix{FStype}(undef,p,p);
        for i=1:m
            ## Compute Fi=f_i(S) in an optimized way
            if (isdiag(S)) # optimize if S is diagonal
                Sd=diag(S);
                local Fid::Vector{FStype}
                if (norm(Sd .- Sd[1])==0) # Optimize further if S is a multiple of identity
                    Fid=fill(ff[i](Sd[1]),p)
                else  # Diagonal but not constant
                    Fid0=Vector{Number}(undef,p)
                    for j=1:p
                        Fid0[j]=ff[i](Sd[j])
                    end
                    Fid=Vector{FStype}(Fid0);
                end
                Fi .= Diagonal(Fid)
            else
                Fi .= ff[i](S);
            end
            mul!(VFi,V,Fi);
            Z[:,:] .+= AA[i]*VFi;
        end
        return Z
    end

""" Internal helper function: returns (TZ,x) where x is a vector of function evaluation of nep.fi[i](λ). TZ is greatest of the types in vector x."""
    function compute_Mder_fi_and_output_type(nep::AbstractSPMF,λ::Number)
        ff=get_fv(nep);
        AA=get_Av(nep)
        x = [f(λ) for f in ff]

        # figure out the return type, as the greatest type of all input
        Tx = promote_typeof(x)
        TA = promote_eltype(AA) # Greatest type of all A-matrices
        TZ = promote_type(TA,Tx)  # output type
        return (TZ,x);
    end

""" Internal helper function: Computes a linear combination with coeffs x of the A-matrices in an
SPMF. The result will be stored in an  AbstractMatrix with eltype T.  """
    function compute_linear_combination_SPMF(nep::AbstractSPMF,::Type{T},x::Vector) where {T<:Number}
        AA=get_Av(nep);
        # Fall back version just sum it up, e.g., for dense matrices
        Z=mapreduce(i->AA[i]*x[i],+,1:length(AA))
        #@assert eltype(Z) <: T
    end
    # Specialized function sparse matrices, where we can exploit sparsity
    function compute_linear_combination_SPMF(nep::SPMF_NEP_sparse,::Type{T},x::Vector) where {T<:Number}
        local Z::SparseMatrixCSC{T,Int}
        if (nep.sparsity_patterns_aligned)
             # Sparsity patterns are aligned, so just change
             # the value entries in the sparse matrix object
             Z = SparseMatrixCSC(nep.A[1].m, nep.A[1].n,
                                 copy(nep.A[1].colptr), copy(nep.A[1].rowval),
                                 copy(convert.(T, nep.A[1].nzval .* x[1])))
             for k = 2:length(nep.A)
                 Z.nzval .+= nep.A[k].nzval .* x[k]
             end
         else
             AA=nep.A;
             # No alignment of sparsity pattern. Naive summing.
             Z=mapreduce(i->AA[i]*x[i],+,1:length(AA))
         end
         return Z
    end

    function compute_Mder(nep::AbstractSPMF,λ::Number)
        TZ,x = compute_Mder_fi_and_output_type(nep,λ)
        # Try to do "smart" summing
        Z = compute_linear_combination_SPMF(nep,TZ,x);
        return Z
     end

    # For higher derivatives
    function compute_Mder(nep::SPMF_NEP{T,Ftype},λ::Number,i::Integer) where {T,Ftype}
        if (i==0)
            return compute_Mder(nep,λ);
        else
            n=size(nep,1);

            # Jordan matrix trick
            k=i+1;
       	    S=diagm(0 => fill(λ,k), -1 => (1:k-1))
            TS=eltype(S);

            # Type logic
            Fλtype=promote_type(TS,Ftype);
            TT=promote_type(Fλtype,eltype(nep.A[1])); # Return type

            # Compute the derivatives of the individual functions with the jordan matrix trick
            fk=Vector{TT}(undef,size(nep.A,1))
            for j=1:size(nep.A,1)
                fk[j] = nep.fi[j](S)[end,1]; # Contains the derivative
            end

            # Try to do "smart" summing
            return compute_linear_combination_SPMF(nep,TT,fk);
        end
    end


    ###########################################################
    # Delay eigenvalue problems - DEP
    #

    """
    type DEP <: AbstractSPMF
    function DEP(AA::Vector{AbstractMatrix} [,tauv::Vector=[0,1.0]])

A `DEP` (Delay Eigenvalue problem) is a problem defined by
the sum
```math
M(λ)=-λI + Σ_i A_i exp(-τ_i λ)
```
where all of the matrices are of size ``n×n``. This type of
NEP describes the stability of time-delay systems.

The construction takes the system matrices ``A_i``, and `tauv` is a vector of the values  ``τ_i``.

# Example:
```julia-repl
julia> A0=randn(3,3); A1=randn(3,3);
julia> tauv=[0,0.2] # Vector with delays
julia> dep=DEP([A0,A1],tauv)
julia> λ=3.0;
julia> M1=compute_Mder(dep,λ)
julia> M2=-λ*I+A0+A1*exp(-tauv[2]*λ)
julia> norm(M1-M2)
0.0
```
"""
    struct DEP{Z<:Real, T<:AbstractMatrix} <: AbstractSPMF{T}
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
                fv[i+1] = S -> one(S)
            else
                # introducing tau is needed to define type stable functions
                tau=nep.tauv[i] #delay of the i-th term
                fv[i+1] = S -> exp(-tau*S)
            end
        end

        return fv
    end

    include("types_poly.jl");

    ###########################################################
    # Rational eigenvalue problem - REP
"""
    function REP(A,roots,poles)

A `REP`-call creates a rational eigenvalue problem. The `REP` is defined by the
sum ``Σ_i A_i s_i(λ)/q_i(λ)``, where i = 0,1,2,..., all of the
matrices are of size ``n×n`` and ``s_i`` and ``q_i`` are polynomials.
The constructor takes the roots and poles as input of polynomials with
normalized highest coefficient. The NEP is defined as
```math
-λI+A_0+A_1\\frac{p(λ)}{q(λ)}
```
where `p` has the roots `roots` and `q` has the roots `poles`.

# Example
```julia-repl
julia> A0=[1 2; 3 4]; A1=[3 4; 5 6];
julia> nep=REP([A0,A1],[1,3], [4,5,6]);
julia> compute_Mder(nep,3)
2×2 Array{Float64,2}:
 Inf  Inf
 Inf  Inf
julia> (λ,x)=quasinewton(nep,v=[1;0])
(-0.3689603779201249 + 0.0im, Complex{Float64}[-2.51824+0.0im, 1.71283+0.0im])
julia> -λ*x+A0*x+A1*x*(λ-1)*(λ-3)/((λ-4)*(λ-5)*(λ-6))
2-element Array{Complex{Float64},1}:
 -2.5055998942313806e-13 + 0.0im
   1.318944953254686e-13 + 0.0im
```
"""
    function REP(AA::Vector,roots::Vector,poles::Vector)

        n=size(AA[1],1)
        A=[one(AA[1]),AA[1],AA[2]];
        roots_float=float.(roots);
        poles_float=float.(poles);
        nep=SPMF_NEP(A,
                     [S-> -S,S->one(S),
                     S-> root_eval(S,poles_float)\root_eval(S,roots_float)]);

        return nep;
    end


    function root_eval(S,a::Vector)
        F=mapreduce(i-> S-a[i]*I, (S1,S2) -> S1*S2, 1:size(a,1))
        return F;
    end




    #######################################################
    ### Represents a projected NEP
"""
    abstract type Proj_NEP <: NEP

`Proj_NEP` represents a projected `NEP`. The projection is defined
as the NEP
```math
N(λ)=W^HM(λ)V
```
where ``M(λ)`` is a base NEP and `W` and `V` rectangular matrices representing
a basis of the projection spaces.
Instances are created with `create_proj_NEP`. See [`create_proj_NEP`](@ref)
for examples.

Any `Proj_NEP` needs to implement two functions to manipulate the projection:

* [`set_projectmatrices!`](@ref): Set matrices `W` and `V`
* [`expand_projectmatrices!`](@ref): Effectively expand the matrices `W` and `V` with one column.

"""
    abstract type Proj_NEP <: NEP end

"""
    pnep=create_proj_NEP(orgnep::ProjectableNEP[,maxsize [,T]])

Create a NEP representing a projected problem ``N(λ)=W^HM(λ)V``,
 where the base NEP is represented by `orgnep`.
The optional parameter `maxsize::Int` determines how large the projected
problem can be and `T` is the Number type used for the projection matrices
(default `ComplexF64`).
These are needed for memory preallocation reasons.
Use [`set_projectmatrices!`](@ref) and [`expand_projectmatrices!`](@ref)
 to specify projection matrices ``V`` and ``W``.

# Example:
The following example illustrates that a projection
of a `NEP` is also a `NEP` and we can for instance
call `compute_Mder` on it:
```julia-repl
julia> nep=nep_gallery("pep0")
julia> V=Matrix(1.0*I,size(nep,1),2);
julia> W=Matrix(1.0*I,size(nep,1),2);
julia> pnep=create_proj_NEP(nep);
julia> set_projectmatrices!(pnep,W,V);
julia> compute_Mder(pnep,3.0)
2×2 Array{Complex{Float64},2}:
  6.08082+0.0im  -5.47481+0.0im
 0.986559+0.0im  -6.98165+0.0im
julia> W'*compute_Mder(nep,3.0)*V  # Gives the same result
2×2 Array{Float64,2}:
 6.08082   -5.47481
 0.986559  -6.98165
```
If you know that you will only use real projection matrices, you
can specify this in at the creation:
```julia-repl
julia> pnep=create_proj_NEP(nep,2,Float64);
julia> set_projectmatrices!(pnep,W,V);
julia> compute_Mder(pnep,3.0)
2×2 Array{Float64,2}:
 6.08082   -5.47481
 0.986559  -6.98165
```
"""
    function create_proj_NEP(orgnep::ProjectableNEP)
         error("Not implemented. All ProjectableNEP have to implement create_proj_NEP.")
    end
    function create_proj_NEP(orgnep::AbstractSPMF,
                             maxsize::Int=min(size(orgnep,1),201),
                             T::Type{<:Number}=ComplexF64)
        return Proj_SPMF_NEP(orgnep,maxsize,T);
    end



"""
    struct Proj_SPMF_NEP <: Proj_NEP

This type represents the (generic) way to project NEPs which are
`AbstractSPMF`. See examples in [`create_proj_NEP`](@ref).
"""
    mutable struct Proj_SPMF_NEP  <: Proj_NEP
        orgnep::AbstractSPMF
        nep_proj::SPMF_NEP; # An instance of the projected NEP
        orgnep_Av::Vector
        orgnep_fv::Vector
        projnep_B_mem::Vector # A vector of matrices
        function Proj_SPMF_NEP(nep::AbstractSPMF,maxsize::Int,
                               subspace_eltype=ComplexF64)
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


            # compute greatest eltype of Av:
            Av_eltype = mapreduce(eltype,promote_type,this.orgnep_Av);

            Bmat_eltype=promote_type(subspace_eltype,Av_eltype)

            # Construct memory for the projected problem matrices
            this.projnep_B_mem=Vector{Matrix{Bmat_eltype}}(undef,size(this.orgnep_fv,1));
            for k=1:size(this.orgnep_fv,1)
                this.projnep_B_mem[k]=zeros(Bmat_eltype,maxsize,maxsize);
            end

            # Reset the projectmatrices
            set_projectmatrices!(this,zeros(size(this.orgnep,1),0),zeros(size(this.orgnep,1),0))



            return this
        end
    end

"""
    set_projectmatrices!(pnep::Proj_NEP,W,V)
Set the projection matrices for the NEP to `W` and `V`, i.e.,
corresponding the NEP: ``N(λ)=W^HM(λ)V``. See also [`create_proj_NEP`](@ref).

# Example
This illustrates if `W` and `V` are vectors of ones, the projected problem
becomes the sum of the rows and columns of the original NEP.
```julia-repl
julia> nep=nep_gallery("pep0")
julia> pnep=create_proj_NEP(nep);
julia> V=ones(200,1);  W=ones(200,1);
julia> set_projectmatrices!(pnep,W,V);
julia> compute_Mder(pnep,0)
1×1 Array{Complex{Float64},2}:
 48.948104019482756 + 0.0im
julia> sum(compute_Mder(nep,0),dims=[1,2])
1×1 Array{Float64,2}:
 48.948104019482955
```

"""
    function set_projectmatrices!(nep::Proj_SPMF_NEP,W,V)
        ## Sets the left and right projected basis and computes
        ## the underlying projected NEP
        m=size(nep.orgnep_Av,1);
        k=size(V,2);
        @assert(k <= size(nep.projnep_B_mem[1],1)) # Don't go outside the prealloc memory
        # For over all i: Compute the expanded matrices
        WT=copy(W')
        B = map(i -> begin
                # Compute the projecte matrix
                nep.projnep_B_mem[i][1:k,1:k]=WT*nep.orgnep_Av[i]*V;
                view(nep.projnep_B_mem[i],1:k,1:k);
                end, 1:m)
        # Keep the sequence of functions for SPMFs
        nep.nep_proj=SPMF_NEP(B,nep.orgnep_fv,check_consistency=false)
    end

"""
    expand_projectmatrices!(nep::Proj_SPMF_NEP, Wnew, Vnew)

The projected NEP is updated by adding the last column of `Wnew` and `Vnew`
to the basis. Note that `Wnew` and `Vnew` contain also the "old" basis vectors.
See also [`create_proj_NEP`](@ref)

# Example:

In the following example you see that the expanded projected problem
has one row and column more, and the leading subblock is the same
as the smaller projected NEP.
```julia-repl
julia> nep=nep_gallery("pep0"); n=size(nep,1);
julia> V=Matrix(1.0*I,n,2); W=Matrix(1.0*I,n,2);
julia> pnep=create_proj_NEP(nep);
julia> set_projectmatrices!(pnep,W,V);
julia> compute_Mder(pnep,0)
2×2 Array{Complex{Float64},2}:
 0.679107+0.0im   -0.50376+0.0im
 0.828413+0.0im  0.0646768+0.0im
julia> Vnew=[V ones(n)]
julia> Wnew=[W ones(n)]
julia> expand_projectmatrices!(pnep,Wnew,Vnew);
julia> compute_Mder(pnep,0)
3×3 Array{Complex{Float64},2}:
 0.679107+0.0im   -0.50376+0.0im  -12.1418+0.0im
 0.828413+0.0im  0.0646768+0.0im   16.3126+0.0im
 -17.1619+0.0im   -10.1628+0.0im   48.9481+0.0im
```

"""
    function expand_projectmatrices!(nep::Proj_SPMF_NEP,Wnew::AbstractMatrix,Vnew::AbstractMatrix)
        w=Wnew[:,end];
        v=Vnew[:,end];
        Av=nep.orgnep_Av;
        k=size(Vnew, 2)-1;
        m = size(nep.orgnep_Av,1);
        @assert(k+1 <= size(nep.projnep_B_mem[1],1)) # Don't go outside the prealloc memory
        WnewT=copy(Wnew[:,1:k]');
        # For over all i: Compute expanded part of matrices:
        B = map(i -> begin
                # Expand the B-matrices
                nep.projnep_B_mem[i][1:k,k+1]=WnewT*Av[i]*v;
                nep.projnep_B_mem[i][k+1,1:(k+1)]=w'*Av[i]*view(Vnew,:,1:(k+1));
                view(nep.projnep_B_mem[i],1:(k+1),1:(k+1));
                end, 1:m)
        # Keep the sequence of functions for SPMFs
        nep.nep_proj=SPMF_NEP(B,nep.orgnep_fv,check_consistency=false)
    end

    # Use delagation to the nep_proj
    compute_MM(nep::Union{Proj_SPMF_NEP},par...)=compute_MM(nep.nep_proj,par...)
    # Use MM to compute Mlincomb for SPMFs
    compute_Mlincomb(nep::Proj_SPMF_NEP,λ::Number,V::AbstractVecOrMat)=
             compute_Mlincomb(nep.nep_proj,λ,V)
    compute_Mlincomb(nep::Proj_SPMF_NEP,λ::Number,V::AbstractVecOrMat,a::Vector)=
             compute_Mlincomb(nep.nep_proj,λ,V,a)
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
    struct GenericSumNEP{NEP1<:NEP,NEP2<:NEP}  <: NEP

See also: [`SumNEP`](@ref), [`SPMFSumNEP`](@ref)
"""
    struct GenericSumNEP{NEP1<:NEP,NEP2<:NEP}  <: NEP
        nep1::NEP1
        nep2::NEP2
    end

"""
    struct SPMFSumNEP{NEP1<:AbstractSPMF,NEP2<:AbstractSPMF}  <: AbstractSPMF{AbstractMatrix}

See also: [`SumNEP`](@ref), [`GenericSumNEP`](@ref)
"""
    struct SPMFSumNEP{NEP1<:AbstractSPMF,NEP2<:AbstractSPMF}  <: AbstractSPMF{AbstractMatrix}
        nep1::NEP1
        nep2::NEP2
    end

    AnySumNEP=Union{GenericSumNEP,SPMFSumNEP};
    # Creator functions for SumNEP
    function SumNEP(nep1::AbstractSPMF,nep2::AbstractSPMF)
        return SPMFSumNEP(nep1,nep2);
    end
    """
    SumNEP{nep1::NEP,nep2::NEP}
    SumNEP{nep1::AbstractSPMF,nep2::AbstractSPMF}

`SumNEP` is a function creating an object that corresponds to a sum of two NEPs,
i.e., if nep is created by `SumNEP` it is defined by
```math
M(λ)=M_1(λ)+M_2(λ)
```
where ``M_1`` and ``M_2`` are defined by `nep1` and `nep2`.

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

See also: [`SPMFSumNEP`](@ref), [`GenericSumNEP`](@ref)

"""
    function SumNEP(nep1::NEP,nep2::NEP)
        return GenericSumNEP(nep1,nep2);
    end

    # Delegate all interface functions
    size(snep::AnySumNEP)=size(snep.nep1)
    size(snep::AnySumNEP,d)=size(snep.nep1,d)
    compute_Mlincomb(nep::AnySumNEP, λ::Number, V::AbstractVecOrMat) =
        (compute_Mlincomb(nep.nep1, λ, V)+compute_Mlincomb(nep.nep2,λ,V))
    compute_Mder(nep::AnySumNEP, λ::Number,i::Int = 0) =
        (compute_Mder(nep.nep1,λ,i)+compute_Mder(nep.nep2,λ,i))
    compute_MM(nep::AnySumNEP, S::AbstractMatrix,V::AbstractMatrix) =
        (compute_MM(nep.nep1,S,V)+compute_MM(nep.nep2,S,V))

    # For SPMFSumNEP, also delegate the get_Av() and get_fv()
    get_Av(nep::SPMFSumNEP) = [get_Av(nep.nep1); get_Av(nep.nep2)]
    get_fv(nep::SPMFSumNEP) = [get_fv(nep.nep1); get_fv(nep.nep2)]


   #######################################################
   ### Functions in common for many NEPs in NEPTypes

   #
"""
    size(nep::Union{DEP,PEP,SPMF_NEP})
    size(nep::Union{DEP,PEP,SPMF_NEP},dim)

Overloads the size functions for NEPs storing size in `nep.n`
"""
    function size(nep::Union{DEP,PEP,SPMF_NEP})
        return (nep.n,nep.n)
    end
    function size(nep::Union{DEP,PEP,SPMF_NEP},dim)
        return nep.n
    end

"""
    issparse(nep)
Returns true/false if the NEP is sparse (if `compute_Mder()` returns sparse)
"""
    function issparse(nep::Union{DEP,PEP,SPMF_NEP})
        return issparse(nep.A[1])
    end


    function get_Av(nep::Union{SPMF_NEP,PEP})
        return nep.A;
    end
    function get_fv(nep::SPMF_NEP)
        return nep.fi;
    end


    include("nep_deflation.jl")
    include("low_rank_nep.jl")
    include("errmeasure.jl")

    # structure exploitation for DEP
    function compute_Mlincomb(nep::DEP,λ::Number,V::AbstractVecOrMat,
                              a::Vector=ones(eltype(V),size(V,2)))
        n=size(V,1); k=size(V,2);
        # Type logic
        TT=promote_type(eltype(V),typeof(λ),eltype(nep.A[1]),eltype(nep.tauv),eltype(a))

        # initialize variables
        z=zeros(TT,n); Vw = Vector{TT}(undef, n); AVw = Vector{TT}(undef, n)
        # with a direct computation one can see that
        # z=-λV[:,1]-V[:,2]+\sum_{j=1}^{length(nep.tauv)} nep.Av[j+1] (V w)
        # where w is the vector with the scaled delays
        for j=1:length(nep.tauv)
            w=Array{TT,1}(exp(-λ*nep.tauv[j])*(-nep.tauv[j]) .^(0:k-1))
            mul!(Vw, V, a.*w)
            mul!(AVw, nep.A[j], Vw)
            z[:] .+= AVw;
        end

        # distinguis the case V is a vector and V is a matrix
        # fix with the proper derivative count
        if (V isa AbstractVector)
            z .-= a[1]*λ*V
        elseif k==1
            z .-= a[1]*λ*V[:]
        else
            z .+= muladd(-λ*a[1],V[:,1],-a[2]*V[:,2])
        end
        return z
    end

    compute_Mlincomb!(nep::DEP,λ::Number,V::AbstractVecOrMat, a::Vector=ones(size(V,2)))=compute_Mlincomb(nep,λ,V, a)

    function compute_Mlincomb!(nep::SPMF_NEP{T,Ftype},
                               λ::Number,
                               V::AbstractVecOrMat,
                               a::Vector=ones(size(V,2))) where {T,Ftype}

        local n,k;
        n=size(V,1);
        k=size(V,2);

    	# we need to assume that the elements of a are different than zero.
    	V[:,findall(x->x==0,a)] .= 0
    	a[findall(x->x==0,a)] .= 1
        local S,TS;
        if (V isa AbstractVector)
            #Vector means just compute matrix vector
            S=λ
            TS = eltype(λ)
        else
            # V matrix means compute linear combination of derivatives. Use
            # scaling trick
            TS= promote_type(typeof(λ),eltype(a));
       	    S=diagm(0 => fill(λ,k), -1 => (a[2:k]./a[1:k-1]).*(1:k-1))
        end

        # Type logic
        Fλtype=promote_type(TS,Ftype);
        TT=promote_type(Fλtype,eltype(V)); # Return type

        z=zeros(TT,n)
        for i=1:size(nep.A,1)
            # Get the function value if V is a vector,
            # otherwise get a vector of scaled derivatives
            Fi1=(V isa AbstractVector) ? nep.fi[i](S) : nep.fi[i](S)[:,1]

            VFi1=V*Fi1
            z .+= nep.A[i]*VFi1
        end

    	return a[1]*reshape(z,n);
    end

    compute_Mlincomb(nep::SPMF_NEP,λ::Number,V::AbstractVecOrMat, a::Vector=ones(size(V,2)))=compute_Mlincomb!(nep,λ,copy(V), copy(a))


    function compute_Mlincomb(
                        nep::PEP,
                        λ::Number,
                        V::AbstractVecOrMat,
                        a::Vector=ones(eltype(V),size(V,2)))

        # Type logic
        TT=promote_type(typeof(λ),eltype(V),eltype(a),eltype(nep.A[1]));

        n=size(nep,1)
        z=zeros(TT,n)
        d=length(nep.A)-1
        k=min(size(V,2),d+1)
        tmp = Vector{TT}(undef, n)

        if iszero(λ)
            for j=0:k-1
                mul!(tmp, nep.A[j+1], view(V,:,j+1))
                z .+= a[j+1] .* factorial(j) .* tmp
            end
        else
            for j=0:k-1
                for i=j:d
                    mul!(tmp, nep.A[i+1], view(V,:,j+1))
                    z .+= a[j+1] .* λ^(i-j) .* (factorial(i)/factorial(i-j)) .* tmp
                end
            end
        end
        return z
    end



    ########## THE FOLLOWING TYPE DerSPMF AND FUNCTION WILL BE DOCUMENTED AND THE CODE BETTER ORGANIZED ##########
    """
        struct DerSPMF{T<:AbstractMatrix,FDtype,TT<:Number} <: AbstractSPMF{T}

    A DerSPMF is a representation of a NEP defined by a Sum of Products of Matrices and Functions [`SPMF_NEP`](@ref). This format makes the execution of [`compute_Mlincomb`](@ref) for the specified `σ` more efficient.
    """
    struct DerSPMF{T<:AbstractMatrix,TT<:Number,FDtype} <: AbstractSPMF{T}
        spmf::AbstractSPMF{T}
        σ::TT
        fD::Matrix{FDtype}
    end

    # implement all the functions
    function size(nep::DerSPMF)
        return size(nep.spmf)
    end
    function size(nep::DerSPMF,dim)
        return size(nep.spmf,dim)
    end
    function get_fv(nep::DerSPMF)
        return get_fv(nep.spmf)
    end
    function get_Av(nep::DerSPMF)
        return get_Av(nep.spmf)
    end

"""
    newspmf=DerSPMF(spmf,σ,m)

Creates a `DerSPMF` representing the NEP `spmf` where the first `m` derivatives of the functions `f_i` in the number
`σ` are precomputed. This will in general speed up the execution of [`compute_Mlincomb`](@ref) for the
selected shift `σ`.

# Parameters

* `spmf` is the original `AbstractSPMF`.

* `σ::Number` specifies where the derivatives will be precomputed.

# Example
```julia-repl
julia> A0=[1 3; 4 5]; A1=[3 4; 5 6];
julia> id_op=S -> one(S) # Note: We use one(S) to be valid both for matrices and scalars
julia> exp_op=S -> exp(S)
julia> nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
julia> m=5 # Precompute 5 derivatives
julia> σ=3.3
julia> dnep=DerSPMF(nep,σ,m)
julia> V=rand(2,m)
julia> z=compute_Mlincomb(dnep,σ,V)  # Same as below ...
2-element Array{Float64,1}:
39.47327945131233
64.31391434063991
julia> z=compute_Mlincomb(nep,σ,V)  # .. but more efficient
2-element Array{Complex{Float64},1}:
 381.008192939787 + 0.0im
 601.379183090249 + 0.0im
```
"""
    function DerSPMF(spmf::AbstractSPMF,σ::Number,m::Int)
          # Compute DD-matrix from get_fv(spmf)

          Av=get_Av(spmf)
          fv=get_fv(spmf)

          TT=promote_type(typeof(σ),eltype(Av[1]))

          # test if the functions fv introduce a super type
          SS=diagm(0=> σ*ones(TT,2m+2),  -1 => (1:2m+1))
          for t=1:length(fv)
              ci=@code_typed(fv[t](SS[1])); TT=promote_type(eltype(ci[end]),TT);
          end

          p=length(fv)
          # matrix for the computation of derivatives

          fD=Matrix{TT}(undef, 2*m+2,p)
          for t=1:p fD[:,t]=fv[t](SS)[:,1] end
          return DerSPMF(spmf,σ,fD);
    end

    function compute_Mlincomb(
                        nep::DerSPMF{T,FDtype},
                        λ::Number,
                        V::AbstractVecOrMat,
                        a::Vector=ones(size(V,2))) where {T,FDtype}

        if λ!=nep.σ
            return compute_Mlincomb(nep.spmf,λ,V,a)
        else

            local n,k,p
            p=size(nep.fD,2)

            if (V isa AbstractVector)
                k=1; n=length(V)
            else
                n,k=size(V)
            end


            Av=get_Av(nep)
            # Type logic
            TT=promote_type(eltype(V),typeof(λ),eltype(Av[1]),eltype(a))
            z=zeros(TT,n)
            VafD=V*(a.*view(nep.fD,1:k,:));
            @inbounds for j=1:p
                z .+= Av[j]*(view(VafD,:,j))
            end
            return z
        end
    end


    include("nep_type_helpers.jl");

end  # End Module
