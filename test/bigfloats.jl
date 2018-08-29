# Unit  tests for bigfloats. Type stability of methods.

if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery

    using Base.Test
end


println("Test NEP-test with BigFloat")
# Create abigfloat NEP
T=BigFloat
A0=ones(T,4,4)-eye(T,4,4)
u=Array{T,1}(1:4); v=u-2;
A1=u*v';
A2=eye(T,4,4);
A2[2,1]=T(pi);

nep=PEP([A0,A1,A2]);

# Run augnewton once to get an accurate solution.
v0=ones(T,4);
λ0=T(0.2)
itercount_star=0;
λiterates_star=zeros(T,100)
myerrmeasure= (λ,v) -> begin
    global itercount_star,λiterates_star
    itercount_star=itercount_star+1;
    λiterates_star[itercount_star]=λ;
    return norm(compute_Mlincomb(nep,λ,v))
end
println("Bigfloat precomputation");
λstar,vstar=augnewton(T, nep,v=v0,λ=λ0,
                      tol=eps(T)*100,
                      errmeasure=myerrmeasure)

bigfloattest=@testset "BigFloat comparison w $T" for T in
    (Float16,Complex32,ComplexF64)

    nep1=PEP(Array{Array{T,2},1}(nep.A))
    global itercount=0;
    global λiterates=zeros(T,100)
    myerrmeasure= (λ,v) -> begin
        global itercount,λiterates
        itercount=itercount+1;
        λiterates[itercount]=λ;
        return abs(T(λstar-λ))

    end


    v0=ones(T,4);
    λ0=T(0.2)
    λ,v=augnewton(T, nep1,v=v0,λ=λ0,
                  tol=eps(real(T))*2,
#                  displaylevel=1,
                  errmeasure=myerrmeasure)
    @test isa(λ,T) # Check output type

    @test abs(T(λ-λstar))<eps(real(T))*8 # Check that we have a solution

    # Check that the kth iterate is the same as the bigfloat iteration
    k=4
    @test abs(T(λiterates_star[k]-λiterates[k]))<eps(real(T))*10

end
