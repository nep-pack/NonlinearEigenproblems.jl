using LinearAlgebra, BenchmarkTools



n=100^2;
m=50;

function mul1!(W,h) W *= h end
function mul2!(W,h) lmul!(h,view(W,:,:)) end
function mul3!(W,h) rmul!(view(W3,:,:),h) end

function mul4!(W,h)
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds W[i,j] = h*W[i,j]
        end
    end
end

V=rand(n,m)
h=rand()

W1=copy(V)
W2=copy(V)
W3=copy(V)
W4=copy(V)


@time begin mul1!(view(W1,:,:),h) end
@time begin mul2!(W2,h) end
@time begin mul3!(W3,h) end
@time begin mul4!(W4,h) end

println(norm(W1-h*V))
println(norm(W2-h*V))
println(norm(W3-h*V))
println(norm(W4-h*V))
