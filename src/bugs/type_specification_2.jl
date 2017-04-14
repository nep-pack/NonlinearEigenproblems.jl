# Show the importance of correct typing from the start


workspace()

# Pkg.add("BenchmarkTools")
using BenchmarkTools;

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

println("Running slow...");
b1=@benchmark slow(n)
show(STDOUT, "text/plain", b1) # Pretty print it

println("\n\nRunning fast...");
b2=@benchmark fast(n)
show(STDOUT, "text/plain", b2) # Pretty print it


println("\n\nCore of the problem is type-casting ...");

z=0;
println("z=0; typeof(z)=",typeof(z));
z+=1+1im;
println("z+=1+1im; typeof(z)=",typeof(z));
