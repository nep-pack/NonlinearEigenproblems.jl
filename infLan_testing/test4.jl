using BenchmarkTools, LinearAlgebra

m=400
n=100000

V=rand(n,m)
h=rand()


function scal_mul!(V,h)
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds V[i,j] = h*V[i,j]
        end
    end
end

function compare(V)
    @btime begin
        V ./= h
    end

    @btime begin
        scal_mul!(V,1/h)
    end

end

compare(V)
