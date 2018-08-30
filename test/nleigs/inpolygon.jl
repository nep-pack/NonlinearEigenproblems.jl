using Test

include(normpath(string(@__DIR__), "..", "..", "src", "nleigs", "inpolygon.jl"))

points = [x + y*im for x=-1:11 for y=-1:11]
polyx = [0, 0, 5, 10, 10]
polyy = [0, 10, 5, 10, 0]

point_inside(p) = inpolygon(real(p), imag(p), polyx, polyy)

function test_inside_count()
    inside = map(point_inside, points)
    @test length(inside) == 13*13
    @test sum(inside) == 96
end

@testset "inpolygon" begin
    # clock-wise
    test_inside_count()

    # counter-clock-wise
    reverse!(polyx)
    reverse!(polyy)
    test_inside_count()

    # edge cases
    @test !point_inside(NaN)
    @test !point_inside(NaN*im)
    @test !point_inside(Inf)
    @test !point_inside(Inf*im)
end
