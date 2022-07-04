#using NonlinearEigenproblemsTest
using Test
using NonlinearEigenproblems.RKHelper

@testset "discretizepolygon" begin
@bench @testset "concave polygon" begin
    poly = [0.0; 0+10im; 5+5im; 10+10im; 10+0im]

    expected_boundary_points = [
        0 + 0im, 0 + 2.2071068im, 0 + 4.4142136im,
        0 + 6.6213203im, 0 + 8.8284271im, 0.73223305 + 9.267767im,
        2.2928932 + 7.7071068im, 3.8535534 + 6.1464466im, 5.4142136 + 5.4142136im,
        6.9748737 + 6.9748737im, 8.5355339 + 8.5355339im, 10 + 9.863961im,
        10 + 7.6568542im, 10 + 5.4497475im, 10 + 3.2426407im,
        10 + 1.0355339im, 8.8284271 + 0im, 6.6213203 + 0im,
        4.4142136 + 0im, 2.2071068 + 0im]

    nr_boundary = 20
    nr_interior = 100
    boundary, interior = discretizepolygon(poly, true, nr_boundary, nr_interior)

    # returned result contains boundary points + polygon vertices + first vertex duplicated
    @test length(boundary) == nr_boundary + length(poly) + 1
    @test boundary[1:nr_boundary] ≈ expected_boundary_points
    @test boundary[nr_boundary+1:nr_boundary+length(poly)] == poly
    @test boundary[nr_boundary+length(poly)+1] == poly[1]

    # points on polygon's boundary should count as inside; expand polygon a tiny bit to counteract rounding errors
    t = eps()*5
    expanded_poly = poly + [-t-t*im, -t+t*im, 0+t*im, t+t*im, t-t*im]
    @test all(p -> inpolygon(real(p), imag(p), real.(expanded_poly), imag.(expanded_poly)), boundary)

    @test length(interior) >= nr_interior
    @test all(p -> inpolygon(real(p), imag(p), real.(poly), imag.(poly)), interior)
end

@bench @testset "narrow Σ" begin
    boundary, interior = discretizepolygon([-10.0-2im, 10-2im, 10+2im, -10+2im], true, 100, 5)
    @test length(boundary) == 100 + 5
    @test length(interior) >= 5
end

@bench @testset "too narrow Σ" begin
    @test_throws ErrorException discretizepolygon([-10.0-0.2im, 10-0.2im, 10+0.2im, -10+0.2im], true, 100, 5)
end

@bench @testset "unit disk" begin
    boundary, interior = discretizepolygon(Vector{ComplexF64}(), true, 100, 100)

    @test length(boundary) == 101
    @test all(isapprox.(abs.(boundary[1:100]), 1.0, rtol = 100*eps()))

    @test length(interior) >= 100
    @test all(abs.(interior) .< 1)
end

@bench @testset "Chebyshev points" begin
    p1 = -2.0 - 1im
    p2 = 2.0 + 1im
    boundary, interior = discretizepolygon([p1, p2], true, 100, 100)

    @test length(boundary) == 102
    @test all(imag((boundary .- p1) / (p2-p1)) .== 0)

    @test length(interior) >= 100
    @test all(imag((interior .- p1) / (p2-p1)) .== 0)
end
end

# To visualize, use something like this:
#using Plots
#import GR
#scatter(real(zz), imag(zz))
#scatter(real(Z), imag(Z))
