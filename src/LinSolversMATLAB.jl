# Linear Solvers using MATLAB as backend

module LinSolversMATLAB
    using MATLAB
    import LinSolvers

    import LinSolvers.EigSolver
    
    export MatlabEigSSolver

    import LinSolvers.eig_solve
    export eig_solve
    
    """
    A linear solver that will call MATLAB eigs()
"""
    type MatlabEigSSolver <: EigSolver
        A
        B

        function MatlabEigSSolver(A,B=spzeros(eltype(A),0))
            this = new(A,B)
            return this
        end
    end


    function eig_solve(solver::MatlabEigSSolver;nev=6,target=0)
        #TODO: The real/complex partition is because of MATLAB.jl limitation in sparse-complex matrices
        eltype_A = eltype(solver.A)
        if !( eltype_A <: Union{Float64, ComplexF64} )
            error("This implementation only supports matrices of type 'Float64' and 'ComplexF64', you have supplied a matrix of type '", eltype_A,"'.")
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

        mat"""
            aa = double($aa_real) + 1i*double($aa_complex);
            bb = double($bb_real) + 1i*double($bb_complex);

            nn = double($n);
            tt = double($t);
            
            [$V,$D] = eigs(aa,bb,nn,tt);
         """

        if nev > 1
            D = diag(D)
        end

        return D,V;

    end



end

