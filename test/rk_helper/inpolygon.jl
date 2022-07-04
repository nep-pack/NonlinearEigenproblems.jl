#using NonlinearEigenproblemsTest
using Test
using NonlinearEigenproblems.RKHelper

@bench @testset "inpolygon" begin
    points = [x + y*im for x=-1:11 for y=-1:11]
    polyx = [0, 0, 5, 10, 10]
    polyy = [0, 10, 5, 10, 0]

    point_inside(p) = inpolygon(real(p), imag(p), polyx, polyy)

    function test_inside_count()
        inside = map(point_inside, points)
        @test length(inside) == 13*13
        @test sum(inside) == 96
    end

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
