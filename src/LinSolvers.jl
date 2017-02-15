module LinSolvers
    export LinSolver
    export DefaultLinSolver
    export BackslashLinSolver    
    export lin_solve
    
    abstract LinSolver;

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
    
    

end
