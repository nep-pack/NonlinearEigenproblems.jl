# Unit test for the lin-solve mehtods (in src/LinSolver.jl)
# The "Gun" problem form gun_native.jl

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems
using Test
using LinearAlgebra

@bench @testset "linsolvers" begin
    TOL = 1e-10;
    nep = nep_gallery("nlevp_native_gun")

    n = size(nep,1);
    λ = 250^2+1im;

    x=ones(n);

    linsolver1=backslash_linsolvercreator(nep,λ);
    y1=lin_solve(linsolver1,x)

    linsolver2=default_linsolvercreator(nep,λ);
    y2=lin_solve(linsolver1,x)

    @test y1 ≈ y2

    y3=compute_Mder(nep,λ)\x;
    @test y1 ≈ y3
    @test y2 ≈ y3
    F=factorize(compute_Mder(nep,λ));
    y4=F\x;
    @test y1 ≈ y4
    @test y2 ≈ y4
    @test y3 ≈ y4

    normy1=norm(y1);
    normy2=norm(y2);
    normy3=norm(y3);
    normy4=norm(y4);

    println("norm(y1,1)=",normy1);
    # To display further digits
    println("normy1-round(128*normy1)/128=",normy1-round(128*normy1)/128, " (for further digits)");
    println("norm(y2,1)=",norm(y2));
    println("normy2-round(128*normy2)/128=",normy2-round(128*normy2)/128);
    println("norm(y3,1)=",norm(y3));
    println("normy3-round(128*normy3)/128=",normy3-round(128*normy3)/128);
    println("norm(y4,1)=",norm(y4));
    println("normy4-round(128*normy4)/128=",normy4-round(128*normy4)/128);

    M=compute_Mder(nep,λ);
    M1norm=norm(M,1);
    r1=M*y1-x;
    r2=M*y2-x;
    r3=M*y3-x;
    r4=M*y4-x;
    @test norm(r1)/M1norm < eps()


    # By first multiplying with M
    z=M*x;
    normz=norm(z)
    println("z=M*x; norm(z)=",normz, " x'*z=",x'*z)
    println("normz-floor(normz*10)/10=",normz-Int(floor(normz*10))/10);
    t=M\z;
    println("norm(M\\z-x)=",norm(t-x));
    @test norm(t-x)/M1norm < eps()

end
