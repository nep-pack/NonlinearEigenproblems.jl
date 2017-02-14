"""
   Linear system solvers
"""
module LinSolvers
    abstract type LinSolver;


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

    function lin_solve(solver::DefaultLinSolver,x::Array)
        return solver.Afact\x
    end
end
