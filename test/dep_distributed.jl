# Run tests for the dep_distributed example

using NonlinearEigenproblems
using Test
using LinearAlgebra

dep_distributed_exact_eigvals = [
   -0.400236388049641 + 0.970633098237807im
   -0.400236388049641 - 0.970633098237807im
   2.726146249832675 + 0.000000000000000im
   -1.955643591177653 + 3.364550574688863im
   -1.955643591177653 - 3.364550574688863im
   4.493937056300693 + 0.000000000000000im
   -1.631513006819252 + 4.555484848248613im
   -1.631513006819252 - 4.555484848248613im
   -1.677320660400946 + 7.496870451838560im
   -1.677320660400946 - 7.496870451838560im]

myerrmeasure=(λ,v) -> minimum(abs.(λ .- dep_distributed_exact_eigvals) ./ abs(λ))

@testset "gallery(dep_distributed" begin
    dep=nep_gallery("dep_distributed");

    n=size(dep,1);

    M1=compute_Mder(dep, dep_distributed_exact_eigvals[1]);
    @test minimum(svdvals(M1))<eps()*100;


    x=Array{Float64,1}(1:n)
    λ=1+1.0im;
    #
    @test norm(compute_Mder(dep,λ)*x-compute_Mlincomb(dep,λ,x))<eps()*10

    # Check derivative computation with finite difference
    zder=compute_Mlincomb(dep,λ,x,[1],1)
    ee=1e-4;
    zder_approx=(compute_Mlincomb(dep,λ+ee,x)-compute_Mlincomb(dep,λ-ee,x))/(2ee)
    @test norm(zder-zder_approx)<ee^2*10




    (λ,V)=iar(dep,σ=3,neigs=5,errmeasure=myerrmeasure, v=ones(n),
              logger=displaylevel,maxit=100,tol=eps()*100)

    @info λ

    @testset "IAR eigval[$i]" for i in 1:length(λ)
        @test estimate_error(myerrmeasure,λ[i],V[:,i])<1e-10
    end

    @bench @testset "Quasinewton eigval[$i]" for i in 1:length(dep_distributed_exact_eigvals[1:3])
        λ0=round(dep_distributed_exact_eigvals[i]*10)/10
        λ,v=quasinewton(ComplexF64,dep,v=ones(n),λ=λ0,
                        #displaylevel=displaylevel,
                        armijo_factor=0.5,maxit=200,
                        errmeasure=myerrmeasure,
                        tol=eps()*100)
        @test abs((dep_distributed_exact_eigvals[i]-λ)/λ)<eps()*100
    end

    @testset "Errors thrown" begin
        @test_throws MethodError nep_gallery("dep_distributed", 15)
        @test_throws MethodError nep_gallery("dep_distributed", t=15)
    end
end
