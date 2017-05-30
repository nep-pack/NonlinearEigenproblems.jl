module LinSolvers
    using MATLAB
    using NEPCore
    using IterativeSolvers
    
    # Linear system of equation solvers
    export LinSolver
    export DefaultLinSolver
    export BackslashLinSolver
    export GMRESLinSolver
    export lin_solve

    export default_linsolvercreator
    export backslash_linsolvercreator
    export gmres_linsolvercreator

    # Eigenvalue solvers
    export EigSolver
    export JuliaEigSolver
    export JuliaEigSSolver
    export MatlabEigSSolver
    export DefaultEigSolver
    export eig_solve

    abstract LinSolver;
    abstract EigSolver;

    import Base.eltype
    export eltype
    import Base.size
    export size
    import Base.*
    export *


##############################################################################
    """
        The linear solver associated with julia factorize()
    """
    type DefaultLinSolver{T_num<:Number, T_mat<:AbstractMatrix, T_fact_mat<:Factorization
} <: LinSolver
        A::T_mat
        Afact::T_fact_mat

        DefaultLinSolver(A::AbstractMatrix{T_num},Afact::Factorization{T_num}) = new(A, Afact)
    end

    function DefaultLinSolver(nep::NEP,λ)
        A=compute_Mder(nep,λ)
        Afact=factorize(A)
        return DefaultLinSolver{eltype(A),typeof(A),typeof(Afact)}(A, Afact)
    end

    function lin_solve{T_num}(solver::DefaultLinSolver{T_num}, x::Array; tol=eps(real(T_num)))
        return solver.Afact\x
    end

"""
   default_linsolvercreator(nep::NEP, λ)\n
      Creates a linear solver of type 'DefaultLinSolver'.
"""
    function default_linsolvercreator(nep::NEP, λ)
        return DefaultLinSolver(nep, λ)
    end


"""
      A linear solver which calls backslash directly (no pre-factorization)
"""
    type BackslashLinSolver{T_num<:Number, T_mat<:AbstractMatrix} <: LinSolver
        A::T_mat
        BackslashLinSolver(A::AbstractMatrix{T_num}) = new(A)
    end
    function BackslashLinSolver(nep::NEP,λ)
        A=compute_Mder(nep,λ)
        return BackslashLinSolver{eltype(A), typeof(A)}(A)
    end

    function lin_solve{T_num}(solver::BackslashLinSolver{T_num}, x::Array; tol=eps(real(T_num)))
        return solver.A\x
    end

"""
   backslash_linsolvercreator(nep::NEP, λ)\n
      Creates a linear solver of type 'BackslashLinSolver'.
"""
    function backslash_linsolvercreator(nep::NEP, λ)
        return BackslashLinSolver(nep, λ)
    end


##############################################################################
"""
      A wrapper that, given a NEP, computes the matrix-vector-product
      using the compute_Mlincomb(...) function for NEPs.
"""
    type Mlincomb_matvec{T_num<:Number, T_nep<:NEP}
        nep::T_nep
        λ::T_num

        Mlincomb_matvec{T_num, T_nep}(nep::T_nep, λ::T_num) = new{T_num,T_nep}(nep, λ)
    end

    function *(M::Mlincomb_matvec, v::AbstractVector)
        return compute_Mlincomb(M.nep, M.λ, v, a=[1])
    end

    function size(M::Mlincomb_matvec, dim=-1)
        return size(M.nep, dim)
    end

    function eltype(M::Mlincomb_matvec)
        return eltype(M.λ)
    end



"""
      A linear solver based on GMRES (built into Julia)
"""
    type GMRESLinSolver{T_num<:Number, T_nep<:NEP} <: LinSolver
        A::Mlincomb_matvec{T_num, T_nep}
        kwargs

        function GMRESLinSolver{T_num, T_nep}(nep::T_nep, λ::T_num, kwargs)
            A = Mlincomb_matvec{T_num, T_nep}(nep, λ)
            new{T_num, T_nep}(A, kwargs)
        end
    end

    function lin_solve{T_num}(solver::GMRESLinSolver{T_num}, x::Array; tol=eps(real(T_num)))
        return gmres(solver.A, x; tol=tol, solver.kwargs...)
    end

"""
   gmres_linsolvercreator(nep::NEP, λ, kwargs=())\n
      Creates a linear solver of type 'GMRESLinSolver'.
      Accepts kwargs which are passed to Julia-built-in-GMRES.
"""
    function gmres_linsolvercreator(nep::NEP, λ, kwargs=())
        return GMRESLinSolver{typeof(λ), typeof(nep)}(nep, λ, kwargs)
    end


##############################################################################
"""
    A linear EP solver that calls Julia's in-built eig() 
"""
    type JuliaEigSolver <: EigSolver
        A
        B

        function JuliaEigSolver(A,B=zeros(eltype(A),0))
            this = new()
            this.A = A
            this.B = B
        
            return this
        end
    end

    function eig_solve(solver::JuliaEigSolver;nev = 1, target = 0)
        if(solver.B != zeros(eltype(solver.A),0))
            D,V = eig(solver.A,solver.B);
        else
            D,V = eig(solver.A);
        end

        #Sort the eigenvalues wrt distance from target, and permute
        I = sortperm(abs(target*ones(size(D,1))-D));
        D = D[I];V = V[:,I];

        #Return the nev closest values to target
        if(nev == 1)
            return D[1],V[:,1]
        end
        
        return D[1:nev],V[:,1:nev];
    end

"""
    A linear EP solve that calls Julia's in-built eigs()
"""
    type JuliaEigSSolver <: EigSolver
        A
        B

        function JuliaEigSSolver(A,B=spzeros(eltype(A),0))
            this = new()
            this.A = A
            this.B = B

            return this
        end

    end

    function eig_solve(solver::JuliaEigSSolver;nev=6,target=0)
        if(solver.B != spzeros(eltype(solver.A),0))
            warn("Julia's eigs() could return erroneous results for GEPs")
            D,V = eigs(solver.A,solver.B;nev=nev,sigma=target)
        else
            D,V = eigs(solver.A;nev=nev,sigma=target)
        end

        if(nev == 1)
            return D[1],V[:,1]
        end

        return D,V
    end


"""
    A linear solver that will call MATLAB eigs()
"""
    type MatlabEigSSolver <: EigSolver
        A
        B

        function MatlabEigSSolver(A,B=spzeros(eltype(A),0))
            this = new()
            this.A = A
            this.B = B

            return this
        end
    end

    
    function eig_solve(solver::MatlabEigSSolver;nev=6,target=0)
        #TODO: The real/complex partition is because of MATLAB.jl limitation in sparse-complex matrices
        eltype_A = eltype(solver.A)
        if !( eltype_A <: Union{Float64, Complex128} )
            error("This implementation only supports matrices of type 'Float64' and 'Complex128', you have supplied a matrix of type '", eltype_A,"'.")
        end
        aa_real = mxarray(real(solver.A))
        aa_complex = mxarray(imag(solver.A))

        if(solver.B == zeros(eltype(solver.A),0))
            solver.B = speye(eltype(solver.A),size(solver.A,1))
        end
        bb_real = mxarray(real(solver.B))
        bb_complex = mxarray(imag(solver.B))

        n = mxarray(nev)
        t = mxarray(target)

        @mput aa_real aa_complex bb_real bb_complex n t
        @matlab begin
            n = double(n);
            aa = double(aa_real) + 1i*double(aa_complex);
            bb = double(bb_real) + 1i*double(bb_complex);

            t = double(t);
            
            (V,D) = eigs(aa,bb,n,t); 
        end

        @mget D V

        if nev > 1
            D = diag(D)
        end

        return D,V;

    end


"""
    Default linear EP solver which calls checks for sparsity and accordingly assigns an appropriate solver
""" 
    type DefaultEigSolver <: EigSolver
        subsolver::EigSolver

        function DefaultEigSolver(A,B=zeros(eltype(A),0))
            this = new()

            if(issparse(A))
                if(B == zeros(eltype(A),0))
                    subsolver = JuliaEigSSolver(A,B);
                else
                    subsolver = MatlabEigSSolver(A,B);
                end
            else
                subsolver = JuliaEigSolver(A,B);
            end
        end
    end

    # Additional constructor that allows for calling the default linear eigensolver with a NEP rather than an explicit matrix
    function DefaultEigSolver(nep::NEP,λ::Number)
        return DefaultEigSolver(compute_Mder(nep,λ))
    end

    function eig_solve(solver::DefaultEigSolver;nev=size(solver.subsolver.A,1),target=0)
        
        return eig_solve(solver.subsolver,nev=nev,target=target)

    end
end
