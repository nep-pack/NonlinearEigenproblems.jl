# Run tests for the dep_distributed example

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/gallery_extra")
push!(LOAD_PATH, pwd()*"/src/gallery_extra/waveguide")

using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers
using IterativeSolvers


using Base.Test



dep=nep_gallery("dep0");
n=size(dep,1);

IAR=@testset "IAR" begin
    @testset "accuracy eigenpairs" begin

        (λ,Q)=iar(dep,σ=3,Neig=5,v=ones(n),
            displaylevel=0,maxit=100,tol=eps()*100)

        @testset "IAR eigval[$i]" for i in 1:length(λ)
            @test norm(compute_Mlincomb(dep,λ[i],Q[:,i]))<eps()*100
        end
    end

    # NOW TEST DIFFERENT ORTHOGONALIZATION METHODS

    # define "by hand" the double GS orthogonalization
    function doubleGS(VV, vv, h)
        h[:]=VV'*vv;    vv[:]=vv-VV*h;
        g=VV'*vv;       vv[:]=vv-VV*g;
        h[] = h[]+g[];  β=norm(vv);
        vv[:]=vv/β
        return β
    end

    orthmethod1 = (VV, vv, h) -> orthogonalize_and_normalize!(VV, vv, h, DGKS )
    orthmethod2 = (VV, vv, h) -> doubleGS(VV, vv, h )
    orthmethod3 = (VV, vv, h) -> orthogonalize_and_normalize!(VV, vv, h, ModifiedGramSchmidt )
    orthmethod4 = (VV, vv, h) -> orthogonalize_and_normalize!(VV, vv, h, ClassicalGramSchmidt )



    @testset "orthogonalization" begin
        # testing different orthogonalization methods
        @testset "DGKS" begin
            (λ,Q,err,V)=iar(dep,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100,orthmethod=orthmethod1)
            @test norm(V'*V-eye(size(V,2)))<1e-6
         end

         @testset "hand_written_doubleGS" begin
             (λ,Q,err,V)=iar(dep,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100,orthmethod=orthmethod2)
             @test norm(V'*V-eye(size(V,2)))<1e-6
          end

          @testset "ModifiedGramSchmidt" begin
              (λ,Q,err,V)=iar(dep,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100,orthmethod=orthmethod3)
              @test norm(V'*V-eye(size(V,2)))<1e-6
           end

           @testset "ClassicalGramSchmidt" begin
               (λ,Q,err,V)=iar(dep,σ=3,Neig=5,v=ones(n),displaylevel=0,maxit=100,tol=eps()*100,orthmethod=orthmethod4)
               @test norm(V'*V-eye(size(V,2)))<1e-6
            end

    end

end
