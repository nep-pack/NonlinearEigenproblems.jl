# This implements squrtm using Schur method.
# Show the importance of correct typing from the start


workspace()

n = 300

function slow(n::Integer)
    for j = 1:n
        for i = 1:n
            temp = 0 #HERE IS THE ISSUE
            for k = j:i
                temp += (1 + 1im)
            end
        end
    end
end

function fast(n::Integer)
    U = zeros(Complex128,n,n);
    for j = 1:n
        for i = 1:n
            temp = zero(Complex128) #HERE IT IS CORRECTLY TYPED
            for k = j:i
                temp += (1 + 1im)
            end
        end
    end
end


@time slow(n)

@time fast(n)

print()

