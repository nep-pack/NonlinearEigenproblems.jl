module LinSolvers
    using ..NEPCore
    using LinearAlgebra
    using SparseArrays
    using SuiteSparse
    using IterativeSolvers
    using LinearMaps
    using ArnoldiMethod

    using ..NEPTypes:DeflatedNEP
    using ..NEPTypes:deflated_nep_compute_Q

    # Linear system of equation solvers
    export LinSolver
    export FactorizeLinSolver
    export BackslashLinSolver
    export GMRESLinSolver
    export DeflatedNEPLinSolver
    export lin_solve

    # Eigenvalue solvers
    export EigSolver
    export DefaultEigSolver
    export EigenEigSolver
    export ArnoldiEigSolver
    export eig_solve

    import Base.eltype
    export eltype
    import Base.size
    export size
    import Base.*
    export *

##############################################################################
"""
    abstract type LinSolver

Structs inheriting from this type are able to solve linear systems associated with
a NEP, for a specific `λ`-value. The most common are direct solvers such as
[`FactorizeLinSolver`](@ref), [`BackslashLinSolver`](@ref) and iterative solvers
such as [`GMRESLinSolver`](@ref).

The LinSolver objects are usually created by the NEP-algorithms through
creator functions, which are passed as parameters.

# Example

The most common usecase is that you want to pass a `linsolvercreator`-function
as a parameter to the NEP-algorithm.
This example shows how you can use solvers based on backslash or `factorize()`.
In the example, `BackslashLinSolver` does not exploit that the system matrix
remains the same throughout the algorithm and is therefore slower.
```julia-repl
julia> nep=nep_gallery("qdep0");
julia> using BenchmarkTools
julia> v0=ones(size(nep,1));
julia> @btime λ,v=quasinewton(nep,λ=-1,v=v0, linsolvercreator=DefaultLinSolverCreator());
  181.017 ms (3779 allocations: 58.61 MiB)
julia> @btime λ,v=quasinewton(nep,λ=-1,v=v0, linsolvercreator=BackslashLinSolverCreator());
  2.040 s (4510 allocations: 553.24 MiB)
```

# Example

The `LinSolver`s are constructed for extendability. This example creates our
own `LinSolver` which uses an explicit formula for the inverse
if the NEP has dimension 2x2.

Create the types and a creator.
```julia
julia> using LinearAlgebra
julia> struct MyLinSolver <: LinSolver
   M::Matrix{ComplexF64}
end
julia> function my_linsolvercreator(nep,λ)
   M=compute_Mder(nep,λ);
   return MyLinSolver(M);
end
```
Explicit import `lin_solve` to show how to solve a linear system.
```julia-repl
julia> import NonlinearEigenproblems.LinSolvers.lin_solve;
julia> function lin_solve(solver::MyLinSolver,b::AbstractVecOrMat;tol=0)
   M=solver.M;
   invM=(1/(det(M)))*[M[2,2] -M[1,2];-M[2,1] M[1,1]]
   return invM*b
end
julia> nep=SPMF_NEP([[1.0 3.0; 4.0 5.0], [2.0 1.0; -1 2.0]], [S->S^2,S->exp(S)])
julia> λ,v=quasinewton(nep,λ=-1,v=[1;1],linsolvercreator=my_linsolvercreator);
```

See also: [`lin_solve`](@ref),
[`FactorizeLinSolver`](@ref), [`FactorizeLinSolver`](@ref),
[`DefaultLinSolverCreator`](@ref),
[`BackslashLinSolver`](@ref), [`BackslashLinSolverCreator`](@ref),
[`GMRESLinSolver`](@ref), [`GMRESLinSolverCreator`](@ref)

"""
    abstract type LinSolver end

##############################################################################
"""
    struct FactorizeLinSolver <: LinSolver

This represents the linear solver associated with julia `factorize()`.
See [`LinSolver`](@ref) and [`FactorizeLinSolverCreator`](@ref) for examples.
"""
    struct FactorizeLinSolver{T} <: LinSolver
        Afact::T
        umfpack_refinements::Int
    end

    function FactorizeLinSolver(nep::NEP, λ, umfpack_refinements)
        A=compute_Mder(nep,λ)
        Afact=factorize(A)
        # Only activate umfpack_refiments on julia >= 1.9  cf #265
        if (isdefined(Afact,:control))
            Afact.control[8]=umfpack_refinements # Set the maximum number of refiments
        end
        return FactorizeLinSolver(Afact, umfpack_refinements)
    end

"""
    lin_solve(solver::LinSolver, b::AbstractVecOrMat; tol=0)

This function solves the linear system represented in `solver::LinSolver` with a
right-hand side `b`. The `tol` kwarg is controlling how accurate the linear
system needs to be solved. A NEP-algorithm will call this solver every
time a linear system associated with `M(λ)` needs to be solved.

This function must be overloaded if a user wants to define their own
way of solving linear systems. See [`LinSolver`](@ref) for examples.
"""
    function lin_solve(solver::FactorizeLinSolver, x::Array; tol = 0)
        solver.Afact \ x
    end


"""
    struct BackslashLinSolver <: LinSolver

This represents a linear solver corresponding to the backslash operator (no pre-factorization).

See also: [`LinSolver`](@ref) and [`BackslashLinSolverCreator`](@ref)
"""
    struct BackslashLinSolver{T_mat} <: LinSolver
        A::T_mat
    end
#    BackslashLinSolver(A::AbstractMatrix{T_num}) where T_num = new(A)

    function BackslashLinSolver(nep::NEP,λ)
        A=compute_Mder(nep,λ)
        return BackslashLinSolver{typeof(A)}(A)
    end

    function lin_solve(solver::BackslashLinSolver{T_mat}, x::Array; tol=0) where T_mat
        return solver.A\x
    end



##############################################################################
"""
    struct GMRESLinSolver <: LinSolver

This represents a solver done with the julia GMRES implementation.

See also: [`LinSolver`](@ref), [`GMRESLinSolverCreator`](@ref)
"""
    struct GMRESLinSolver{T_num<:Number, T_kwargs} <: LinSolver
        A::LinearMap{T_num}
        kwargs::T_kwargs

        function GMRESLinSolver{T_num}(nep, λ::T_num, kwargs::T_kwargs) where {T_num<:Number, T_kwargs}
            f = v -> compute_Mlincomb(nep,λ,v)
            A = LinearMap{T_num}(f, size(nep,1), ismutating=false)
            new{T_num, T_kwargs}(A, kwargs)
        end

    end


    function lin_solve(solver::GMRESLinSolver{T_num,T_kwargs}, b::Vector{T_num}; tol=eps(real(T_num))) where {T_num,T_kwargs}
        v = zero(b)
        retval=gmres!(v, solver.A, b; tol=tol, solver.kwargs...)
        return v
    end



##############################################################################
"""
    struct DeflatedNEPLinSolver <: LinSolver
The `DeflatedNEPLinSolver` solves the linear system connected with a deflated NEP.
The extended linear system is solveed with a Schur complement strategy, recycling
the linear solver of the original NEP. such that
```math
[M, U; X^T, 0][v1; v2] = [b1; b2]
```
NB1: The implementation assumes minimality index = 1.
NB2: The Schur complement is explicitly formed. Hence it is only efficient for a
few deflated eigenvalues.
See also: [`LinSolver`](@ref), [`DeflatedNEPLinSolverCreator`](@ref),
[`deflate_eigpair`](@ref)
# References
* C. Effenberger, Robust successive computation of eigenpairs for nonlinear eigenvalue problems. SIAM J. Matrix Anal. Appl. 34, 3 (2013), pp. 1231-1256.
"""
        struct DeflatedNEPLinSolver{T_NEP, T_num, T_LinSolver} <: LinSolver
            deflated_nep::T_NEP
            λ::T_num
            orglinsolver::T_LinSolver
        end


        function DeflatedNEPLinSolver(nep::DeflatedNEP, λ, orglinsolver::LinSolver)
            DeflatedNEPLinSolver{typeof(nep), typeof(λ), typeof(orglinsolver)}(nep, λ, orglinsolver)
        end


        function lin_solve(solver::DeflatedNEPLinSolver, b; tol=0)
            deflated_nep = solver.deflated_nep
            λ = solver.λ
            orglinsolver = solver.orglinsolver
            orgnep = deflated_nep.orgnep

            X = deflated_nep.V0
            Λ = deflated_nep.S0

            n = size(orgnep,1)
            m = size(Λ,1)
            T = eltype(b)
            v = zeros(T,n+m)

            v1 = view(v, 1:n)
            v2 = view(v, (n+1):(n+m))
            b1 = b[1:n]
            b2 = b[(n+1):(n+m)]
            U = deflated_nep_compute_Q(deflated_nep, λ, 0)

            # Precompute some reused entities
            b1tilde = lin_solve(orglinsolver, b1, tol=tol) # b1tilde = M^{-1}b1
            Z::Matrix{T} = zeros(T,n,m)
            for i = 1:m
                Z[:,i] = lin_solve(orglinsolver, vec(U[:,i]), tol=tol) # Z = M^{-1}U
            end
            S = -X'*Z #Schur complement
            v2[:] = S\(b2 - X'*b1tilde)
            v1[:] = b1tilde - Z*v2

            return v
        end



##############################################################################
"""
    abstract type EigSolver

Structs inheriting from this type are able to solve linear eigenvalue problems
arising in certain methods, such as, e.g., `mslp`, `sgiter`,
and `polyeig`.

The `EigSolver` objects are passed as types to the NEP-algorithms,
which uses it to dispatch the correct version of the function [`eig_solve`](@ref).

# Example

The most common usecase is that you do not want to specify anything in particular, since
the [`DefaultEigSolver`](@ref) will use a dense or a sparse method depending on you problem.
However, this example shows how you can force `mslp` to use the sparse solver.
```julia-repl
julia> nep=nep_gallery("qdep0");
julia> λ,v = mslp(nep, eigsolvertype=ArnoldiEigSolver);
julia> norm(compute_Mlincomb(nep,λ,v))
9.323110647141726e-16
```

# Example

The `EigSolver`s are constructed for extendability. As an illustartion this example
creates a naive `EigSolver` which casts the problem to a standard linear eigenproblem
and calls the built-in function to solve it.

Create the types and a creator.
```julia-repl
julia> struct MyEigSolver <: EigSolver
   A
   E
   function MyEigSolver(A,E)
      return new(A,E)
   end
end

julia> import NonlinearEigenproblems.LinSolvers.eig_solve;
julia> function eig_solve(solver::MyEigSolver;nev = 1, target = 0)
   M = solver.E \\ solver.A
   eig = eigen(M)
   i = argmin(abs.(eig.values))
   return eig.values[i], eig.vectors[:,i]
end
julia> nep=nep_gallery("dep0", 50);
julia> λ,v = mslp(nep, eigsolvertype=MyEigSolver, tol=1e-5);
julia> norm(compute_Mlincomb(nep,λ,v))
3.0777795031319117e-10
```

See also: [`eig_solve`](@ref),
[`DefaultEigSolver`](@ref), [`EigenEigSolver`](@ref),
[`ArnoldiEigSolver`](@ref), [`eig_solve`](@ref)

"""
    abstract type EigSolver end


##############################################################################
"""
    struct EigenEigSolver <: EigSolver

A linear eigenvalueproblem solver that calls Julia's in-built eigen()

Constructed as `EigenEigSolver(A, [B,])`, and solves the problem
```math
Ax = λBx
```
The paramter `B` is optional an default is indentity, for which a standard linear
eigenproblem is solved.

See also: [`EigSolver`](@ref) and [`eig_solve`](@ref)
"""
    struct EigenEigSolver{T_A,T_B} <: EigSolver
        A::T_A
        B::T_B

        function EigenEigSolver(A)
            return new{typeof(A),Missing}(A,missing)
        end
        function EigenEigSolver(A,B)
            return new{typeof(A),typeof(B)}(A,B)
        end
    end



    function eig_solve(solver::EigenEigSolver; nev = 1, target = 0)
        D,V = inner_eig_solve(solver)

        #Sort the eigenvalues wrt distance from target, and permute
        I = sortperm(abs.(target*ones(size(D,1))-D));
        D = D[I];V = V[:,I];
        return D[1:nev],V[:,1:nev];
    end

    function inner_eig_solve(solver::EigenEigSolver{T_A,T_B}) where {T_A, T_B}
        D,V = eigen(solver.A,solver.B)
    end
    function inner_eig_solve(solver::EigenEigSolver{T_A,T_B}) where {T_A, T_B<:Missing}
        D,V = eigen(solver.A)
    end

"""
    struct ArnoldiEigSolver <: EigSolver

A linear eigenproblem solver for large and sparse problems that calls
the Arnoldi method implemented in the Julia package ArnoldiMethod.jl.

Constructed as `ArnoldiEigSolver(A, [B,])`, and solves the problem
```math
Ax = λBx
```
The paramter `B` is optional an default is indentity, for which a standard linear
eigenproblem is solved.

See also: [`EigSolver`](@ref) and [`eig_solve`](@ref)
"""
    struct ArnoldiEigSolver{T_A,T_B} <: EigSolver
        A::T_A
        B::T_B

        function ArnoldiEigSolver(A)
            return new{typeof(A),Missing}(A,missing)
        end
        function ArnoldiEigSolver(A,B)
            return new{typeof(A),typeof(B)}(A,B)
        end

    end


    function eig_solve(solver::ArnoldiEigSolver; nev=6, target=0)
        D,V = inner_eigs_solve(solver, nev, target)
        return D,V
    end

    function inner_eigs_solve(solver::ArnoldiEigSolver{T_A,T_B}, nev, target) where {T_A, T_B}
        if (T_B <: Missing)
            C=target*I-solver.A;
            Cfact=factorize(C);
            Atransformed=LinearMap{eltype(Cfact)}(x->Cfact\x,
                                   size(solver.A,1),size(solver.A,1));
        else

            C=target*solver.B-solver.A;
            Cfact=factorize(C);
            Atransformed=LinearMap{eltype(Cfact)}(x->Cfact\(solver.B*x),
                                   size(solver.A,1),size(solver.A,1));
        end
        # Call restarted Arnoldi
        decomp, history = partialschur(Atransformed, nev=nev, tol=1e-10, which=LM());
        D0, V = partialeigen(decomp)
        # Sort. So we don't depend on partialschur to do be sorted already
        IJ=sortperm(-abs.(D0))
        # And reverse transformation
        D = target .- inv.(D0) # Reverse transformation
        return D[IJ[1:nev]],V[:,IJ[1:nev]]
    end




"""
    struct DefaultEigSolver <: EigSolver

A linear eigenvalueproblem solver that calls checks for sparsity and accordingly
assigns an appropriate solver.

See also: [`EigSolver`](@ref), [`eig_solve`](@ref), [`EigenEigSolver`](@ref), [`ArnoldiEigSolver`](@ref)
"""
    struct DefaultEigSolver{T_sub} <: EigSolver
        subsolver::T_sub

        function DefaultEigSolver(A,B)
            local subsolver
            if(issparse(A))
                subsolver = ArnoldiEigSolver(A,B)
            else
                subsolver = EigenEigSolver(A,B)
            end
            return new{typeof(subsolver)}(subsolver)
        end
        function DefaultEigSolver(A)
            local subsolver
            if(issparse(A))
                subsolver = ArnoldiEigSolver(A)
            else
                subsolver = EigenEigSolver(A)
            end
            return new{typeof(subsolver)}(subsolver)
        end
    end


"""
    eig_solve(solver::EigSolver; [nev,] [target,])

This function solves the linear eigenvalue problem represented in `solver::EigSolver`.
The `nev` kwarg is controlling the number of eigenvalues aimed for, and `target`
specifies around which point the eigenvalues are computed. The former has a defalut value
equalt to the seize of the problem, and the latter has a defalut value 0.

Return values are of the form (Vector, Matrix) where the former contains the eigenvalues
and the latter the eigenvectors.

This function must be overloaded if a user wants to define their own
way of solving linear eigenvalue problems. See [`EigSolver`](@ref) for examples.
"""
    function eig_solve(solver::DefaultEigSolver;nev=size(solver.subsolver.A,1),target=0)
        return eig_solve(solver.subsolver,nev=nev,target=target)
    end

    include("LinSolverCreators.jl");

end
