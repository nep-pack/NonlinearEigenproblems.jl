using Base.Test

include(normpath(string(@__DIR__), "..", "..", "src", "nleigs", "discretizepolygon.jl"))
include(normpath(string(@__DIR__), "..", "..", "src", "nleigs", "inpolygon.jl"))

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

@testset "discretizepolygon" begin
    # returned result contains boundary points + polygon vertices + first vertex duplicated
    @test length(boundary) == nr_boundary + length(poly) + 1
    @test boundary[1:nr_boundary] ≈ expected_boundary_points
    @test boundary[nr_boundary+1:nr_boundary+length(poly)] == poly
    @test boundary[nr_boundary+length(poly)+1] == poly[1]

    # points on polygon's boundary should count as inside; expand polygon a tiny bit to counteract rounding errors
    t = eps()*5;
    expanded_poly = poly + [-t-t*im, -t+t*im, 0+t*im, t+t*im, t-t*im]
    @test all(p -> inpolygon(real(p), imag(p), real.(expanded_poly), imag.(expanded_poly)), boundary)

    @test length(interior) >= nr_interior
    @test all(p -> inpolygon(real(p), imag(p), real.(poly), imag.(poly)), interior)
end

# To visualize:
if false
    using Plots
    import GR
    scatter(real(zz), imag(zz))
    scatter(real(Z), imag(Z))
end
