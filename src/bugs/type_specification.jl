# This example is taken from the Julia manual:
# http://docs.julialang.org/en/stable/manual/performance-tips/#type-declarations


workspace()

type MyType{T<:AbstractFloat}
    a::T
end

func(m::MyType) = m.a+1

println("First version, type is concrete")
code_llvm(func,(MyType{Float64},))

println("\n\n\nSecond version, type is abstract")
code_llvm(func,(MyType{AbstractFloat},))

println("\n\n\nThird version, type is not given")
code_llvm(func,(MyType,))

