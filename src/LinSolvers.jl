module LinSolvers
    export LinSolver
    export DefaultLinSolver
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
end
