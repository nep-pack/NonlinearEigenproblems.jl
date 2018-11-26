using BenchmarkTools, LinearAlgebra

m=400
n=10000

V=rand(n,m)
W=rand(n,m)
h=rand()

function mul_and_sub!(V,W,h)

    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds V[i,j] -= h*W[i,j]
        end
    end

end

function compare(V,W,h)
    @btime begin
        V[:,:] .-= h*W;
    end

    @btime begin
        mul_and_sub!(V,W,h)
    end
    VV=copy(V)
    V.-= h*W;
    mul_and_sub!(VV,W,h)
    norm(V-VV)
end

compare(V,W,h)
