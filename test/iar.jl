using NonlinearEigenproblems
using Test
using IterativeSolvers
using LinearAlgebra

# The user can create his own orthogonalization function to use in IAR
function doubleGS_function!(VV, vv, h)
    h[:]=VV'*vv; vv[:]=vv-VV*h; g=VV'*vv; vv[:]=vv-VV*g;
    h[:] = h[:]+g[:]; β=norm(vv); vv[:]=vv/β; return β
end
# Then it is needed to create a type to access to this function
abstract type DoubleGS <: IterativeSolvers.OrthogonalizationMethod end
# And then introduce a function dispatch for this new type in order to use
# the defined orthogonalization function
import IterativeSolvers.orthogonalize_and_normalize!
function orthogonalize_and_normalize!(V,v,h,::DoubleGS)
    doubleGS_function!(V, v, h) end

@testset "IAR" begin
    dep=nep_gallery("dep0");
    n=size(dep,1);

    @bench @testset "accuracy eigenpairs" begin
        (λ,Q)=iar(dep,σ=1.1,neigs=5,v=ones(n),
                  maxit=100,tol=eps()*100,errmeasure=ResidualErrmeasure(dep));
        verify_lambdas(5, dep, λ, Q, eps()*100)
    end

    @bench @testset "Solve by projection" begin
        (λ,Q)=iar(dep,σ=1.1,neigs=5,v=ones(n),
                  maxit=100,tol=eps()*100,errmeasure=ResidualErrmeasure(dep), proj_solve=true);
        verify_lambdas(5, dep, λ, Q, eps()*100)
    end

    @testset "Compute as many eigenpairs as possible (neigs=Inf)" begin
        (λ,Q)=iar(dep,σ=1.1,neigs=Inf,v=ones(n),
                  maxit=38,tol=eps()*100);
        verify_lambdas(6, dep, λ, Q, eps()*200)
    end

    @testset "orthogonalization" begin
        # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS

        @bench @testset "DGKS" begin
            (λ,Q,V)=iar(dep,orthmethod=DGKS(),σ=1.1,neigs=5,v=ones(n),maxit=100,tol=eps()*100)
            @test opnorm(V'*V - I) < 1e-6
        end

        @bench @testset "User provided doubleGS" begin
            (λ,Q,V)=iar(dep,orthmethod=DoubleGS(),σ=1.1,neigs=5,v=ones(n),maxit=100,tol=eps()*100)
            @test opnorm(V'*V - I) < 1e-6
        end

        @bench @testset "ModifiedGramSchmidt" begin
            (λ,Q,V)=iar(dep,orthmethod=ModifiedGramSchmidt(),σ=1.1,neigs=5,v=ones(n),maxit=100,tol=eps()*100)
            @test opnorm(V'*V - I) < 1e-6
        end

        @bench @testset "ClassicalGramSchmidt" begin
            (λ,Q,V)=iar(dep,orthmethod=ClassicalGramSchmidt(),σ=1.1,neigs=5,v=ones(n),maxit=100,tol=eps()*100)
            @test opnorm(V'*V - I) < 1e-6
        end
    end

    @testset "Errors thrown" begin
        np=100;
        dep=nep_gallery("dep0",np);
        @test_throws NEPCore.NoConvergenceException (λ,Q)=iar(dep,σ=1.1,neigs=6,v=ones(np),
                  maxit=7,tol=eps()*100);
    end

end
