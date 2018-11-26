using BenchmarkTools, LinearAlgebra

m=400
n=1000

V=rand(n,m)


function mat_zero!(V)
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds V[i,j] = 0
        end
    end
end

function compare(V)
    @btime begin
        V[:,:] = zero(V);
    end


    @btime begin
        mat_zero!(V)
    end

end

compare(V)
