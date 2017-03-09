module LinSolvers
    using MATLAB
    # Linear system of equation solvers
    export LinSolver
    export DefaultLinSolver
    export BackslashLinSolver    
    export lin_solve

    # Eigenvalue solvers
    export EigSolver
    export JuliaEigSolver
    export JuliaEigSSolver
    export MatlabEigSSolver
    export DefaultEigSolver
    export eig_solve

    abstract LinSolver;
    abstract EigSolver;

##############################################################################
    """
        The linear solver associated with julia factorize()
    """
    type DefaultLinSolver <: LinSolver
        A
        Afact
        function DefaultLinSolver(A)
            this=new()
            this.Afact=factorize(A)
            return this
        end
    end
    function lin_solve(solver::DefaultLinSolver,x::Array;tol=eps())
        return solver.Afact\x
    end
"""
      A linear solver which calls backslash directly (no pre-factorization)
"""
    type BackslashLinSolver <: LinSolver
        A
        function BackslashLinSolver(A)
            this=new(A)
            return this
        end
    end

    function lin_solve(solver::BackslashLinSolver,x::Array;tol=eps())
        return solver.A\x
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

        aa = mxarray(solver.A)

        #is_gep = mxarray(0);
        if(solver.B != zeros(eltype(solver.A),0))
            bb = mxarray(solver.B)
            #is_gep = mxarray(1);
        else 
            solver.B = speye(eltype(solver.A),size(solver.A,1))
            bb = mxarray(solver.B)        
        end

        n = mxarray(nev)
        t = mxarray(target)

        @mput aa bb n t 
        @matlab begin
            n = double(n);
            aa = double(aa);
            bb = double(bb);

            t = double(t);
            #is_gep = double(is_gep);
            
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
        A
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

    function eig_solve(solver::DefaultEigSolver;nev=size(solver.subsolver.A,1),target=0)
        
        return eig_solve(solver.subsolver,nev=nev,target=target)

    end
end
