#LEJABAGBY generate Leja-Bagby points (a,b) on (A,B), with
#  scaling factors beta such that the uniform norm on the control set C is
#  1. Greedy search for a minimum is performed on B. If keepA is true then
#  the points in the output a will be exactly those of A, otherwise the
#  points in a are also chosen via greedy search on A. If forceInf is a
#  positive integer, the first forceInf poles in b will be infinity.
function lejabagby(A::AbstractVector{Complex{T}}, B::AbstractVector{T}, C::AbstractVector{Complex{T}}, m::Int, keepA::Bool=false, forceInf::Int=0) where T<:Real
    if minimum(abs.(B)) < 1e-9
        warn("There is at least one pole candidate in B being nearby zero. Consider shifting your problem for stability.")
    end

    a = [A[1]]
    b = [forceInf > 0 ? Inf : B[1]]
    beta = [1.]

    sA = ones(A)
    sB = ones(eltype(A), size(B))
    sC = ones(eltype(A), size(C))

    for j = 1:m-1
        sA .*= ((A-a[j]) ./ (1 - A/b[j]));
        sB .*= ((B-a[j]) ./ (1 - B/b[j]));
        sC .*= ((C-a[j]) ./ (1 - C/b[j]));

        push!(a, A[keepA ? j+1 : indmax(abs.(sA))])
        push!(b, forceInf > j ? Inf : B[indmin(abs.(sB))])
        push!(beta, maximum(abs.(sC)))

        # treat single point case
        if beta[j+1] < eps()
            beta[j+1] = 1
        end

        sA /= beta[j+1]
        sB /= beta[j+1]
        sC /= beta[j+1]
    end

    return a, b, beta
end
