using NonlinearEigenproblems
using Test
using LinearAlgebra
using IterativeSolvers

# The user can create his own orthogonalization function to use in IAR
function doubleGS_function!(VV, vv, h)
    h[:]=VV'*vv; vv[:]=vv-VV*h; g=VV'*vv; vv[:]=vv-VV*g;
    h[] = h[]+g[]; β=norm(vv); vv[:]=vv/β; return β
end
# Then it is needed to create a type to access to this function
struct DoubleGS <: IterativeSolvers.OrthogonalizationMethod end
# And then introduce a function dispatch for this new type in order to use
# the defined orthogonalization function
import IterativeSolvers.orthogonalize_and_normalize!
function orthogonalize_and_normalize!(V,v,h,::DoubleGS)
    doubleGS_function!(V, v, h) end

@testset "TIAR" begin
    n=100;
    dep=nep_gallery("dep0",n);

    @bench @testset "accuracy eigenpairs" begin
        (λ,Q)=tiar(dep,σ=1.1,γ=3,neigs=2,v=ones(n),maxit=50,tol=eps()*100,errmeasure=ResidualErrmeasure(dep));
        verify_lambdas(2, dep, λ, Q, eps()*100)
    end

    @testset "Compute as many eigenpairs as possible (neigs=Inf)" begin
        (λ,Q)=tiar(dep,σ=1.1,γ=3,neigs=Inf,v=ones(n),maxit=50,tol=eps()*100,errmeasure=ResidualErrmeasure(dep));
        verify_lambdas(7, dep, λ, Q, eps()*100)
    end

    @testset "orthogonalization" begin

    # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS
    @bench @testset "DGKS" begin
        (λ,Q,Z)=tiar(dep,σ=1.1,γ=3,neigs=4,v=ones(n),maxit=50,tol=eps()*100,errmeasure=ResidualErrmeasure(dep),orthmethod=DGKS());
        @test opnorm(Z'*Z - I) < 1e-6
     end

     @bench @testset "User provided doubleGS" begin
         (λ,Q,Z)=tiar(dep,σ=1.1,γ=3,neigs=4,v=ones(n),maxit=50,tol=eps()*100,orthmethod=DoubleGS());
         @test opnorm(Z'*Z - I) < 1e-6
      end

      @bench @testset "ModifiedGramSchmidt" begin
          (λ,Q,Z)=tiar(dep,σ=1.1,γ=3,neigs=4,v=ones(n),maxit=50,tol=eps()*100,orthmethod=ModifiedGramSchmidt());
          @test opnorm(Z'*Z - I) < 1e-6
      end

       @bench @testset "ClassicalGramSchmidt" begin
           (λ,Q,Z)=tiar(dep,σ=1.1,γ=3,neigs=4,v=ones(n),maxit=50,tol=eps()*100,orthmethod=ClassicalGramSchmidt());
           @test opnorm(Z'*Z - I) < 1e-6
       end
    end

    # iar and tiar ara mathematically equivalent it maxit<<nep
    # verify the equivalence numerically
    @bench @testset "equivalance with IAR" begin
        (λ_tiar,Q_tiar)=tiar(dep,σ=1.1,γ=3,neigs=2,v=ones(n),maxit=50,tol=eps()*100);
        (λ_iar,Q_iar)=   iar(dep,σ=1.1,γ=3,neigs=2,v=ones(n),maxit=50,tol=eps()*100);

        # Relaxed to only test abs-values to avoid problems if iar / tiar finds
        # the other conjugate solution
        for k=1:size(λ_tiar,1);
            @test minimum(abs.(abs(λ_iar[k]) .- abs.(λ_tiar))) < 1e-6
            @test minimum(map(j-> norm(abs.(Q_tiar[:,k]/Q_tiar[1,k])-abs.(Q_iar[:,j]/Q_iar[1,k])),1:size(λ_iar,1)))<1e-6
        end

    end
    @bench @testset "Solve by projection" begin
        np=1000;

        depp=nep_gallery("dep0",np);
        nn=opnorm(compute_Mder(depp,0));
        errmeasure= (λ,v) -> norm(compute_Mlincomb(depp,λ,v))/nn;

        λ,Q = tiar(depp, σ=0, γ=3, neigs=3, v=ones(np), maxit=50,
                   tol=sqrt(eps()), check_error_every=3,
                   proj_solve=true, inner_solver_method=IARInnerSolver(),
                   errmeasure=errmeasure);

        @test errmeasure(λ[1],Q[:,1])<sqrt(eps())*10
    end

    @testset "Errors thrown" begin
        np=100;
        dep=nep_gallery("dep0",np);
        @test_throws NEPCore.NoConvergenceException (λ,Q)=tiar(dep,σ=2.0,γ=3,neigs=4,v=ones(np),maxit=5,tol=eps()*100,errmeasure=ResidualErrmeasure(dep));
    end

end
