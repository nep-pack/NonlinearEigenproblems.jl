# Run tests for the dep_distributed example

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using IterativeSolvers
using Base.Test



dep=nep_gallery("dep_distributed");
n=size(dep,1);
exact_eigvals=[
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

@testset "gallery(dep_distributed" begin
    M1=compute_Mder(dep,exact_eigvals[1]);
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


    myerrmeasure=(λ,v) -> begin
        global exact_eigvals,dep;
        return minimum(abs.(λ-exact_eigvals)./abs(λ))
        #return norm(compute_Mlincomb(dep,λ,v))/norm(compute_Mder(dep,λ))
    end

    (λ,V)=iar(dep,σ=3,Neig=5,errmeasure=myerrmeasure, v=ones(n),
              displaylevel=1,maxit=100,tol=eps()*100)
    println(λ)

    @testset "IAR eigval[$i]" for i in 1:length(λ)
        @test myerrmeasure(λ[i],V[:,i])<1e-10
    end


    @testset "Quasinewton eigval[$i]" for i in 1:length(exact_eigvals)
        λ0=round(exact_eigvals[i]*10)/10
        λ,v=quasinewton(Complex128,dep,v=ones(n),λ=λ0,
                        #displaylevel=1,
                        armijo_factor=0.5,maxit=200,
                        errmeasure=myerrmeasure,
                        tol=eps()*100)
        @test abs((exact_eigvals[i]-λ)/λ)<eps()*100
    end




end
