
    export sg_iteration



    #############################################################################
    # Safeguarded iteration
    #############################################################################
    sg_iteration(nep::NEP;params...)=sg_iteration(Complex128,nep;params...)
"""
    Safeguarded iteration
"""
    function sg_iteration{T}(::Type{T}, nep::NEP;
                        errmeasure::Function =
                        default_errmeasure(nep::NEP),
                        tolerance_outer = eps(real(T))*10000,
                        tolerance_inner = tolerance_outer/100,
                        λ = zero(T),
                        v = randn(nep.n),
                        displaylevel=0,
                        maxit = 10,
                        eigsolvertype::DataType=DefaultEigSolver)




    #The main program

    println("Running safeguarded iteration, initial approximation of λ: ",λ)
    λ_m::T = T(λ)
    v_m::Array{T,1} = Array{T,1}(v)
    err = T(Inf)


    for k=1:maxit

        #Newton as an inner loop for the Rayleigh quotient
        println("k: ",k,"λ_m: ",λ_m);
        λ_m = compute_rf(nep,v_m,y=v_m, λ0=λ_m,TOL=tolerance_inner);
        println("k: ",k,"λ_m: ",λ_m);

        # Find closest eigenvalue to λ
            # and the corresponding eigenvector
        solver = eigsolvertype(nep,λ_m);
        d,v_m = eig_solve(solver,target=λ_m,nev=1);

        
        # Stopping criterion using err_measure()
        err=errmeasure(λ_m,v_m);
        if (err < tolerance_outer)
            println("Solution found at step ",k);
            println("λ = ",λ_m);
            return (λ_m,v_m)
        end

        println("Error measure at step ",k," = ", err);
        
        abs(d - λ_m)
        
    end
                
    err=errmeasure(λ_m,v_m)
    msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ_m,v_m,err,msg))

end
