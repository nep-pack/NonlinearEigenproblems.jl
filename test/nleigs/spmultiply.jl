using Base.Test

include(normpath(string(@__DIR__), "..", "..", "src", "nleigs", "spmultiply.jl"))

# Run expression 'ex' 'iterations' times and report median time and memory usage
macro timem(name, iterations, ex)
    quote
        t = [collect((@timed $(esc(ex)))[2:3]) for x=1:$iterations][2:end]
        m = median(reduce(hcat, t), 2)
        @printf("%s%7.3f ms (%.0f kb)\n", $(esc(name)), m[1]*1e3, m[2]/1e3)
    end
end

@testset "spmultiply" begin
    n = 2_000_000
    p = 2e-5
    A = sprand(n, n, p)
    v = sprand(n, 8/n)

    @timem "Regular multiply " 10 A*v
    @timem "spmultiply       " 10 spmultiply(A, v)

    @test norm(A*v - spmultiply(A, v)) == 0
end
