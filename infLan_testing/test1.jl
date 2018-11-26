using BenchmarkTools

m=400
n=100000

V=rand(m,n)
W=rand(m,n)

function compare(V,W)
    @btime begin
        α=sum(V.*W)
    end

    @btime begin
        β=mat_sum(V,W)
    end
    α=sum(V.*W)
    β=mat_sum(V,W)
    abs(α-β)


end

function mat_sum(V,W)
    TT=eltype(V)
    β::TT=0
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds β += V[i,j]*W[i,j]
        end
    end
    return β
end

compare(V,W)
