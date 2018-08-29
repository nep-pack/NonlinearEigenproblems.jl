example = 1

n = 10;
Random.seed!(0)

if example == 1
println("Example 1, force cast to sparse - Not working")
  A = rand(n,n)
  A = sparse(A)

  b = rand(n)
  b = sparse(b)
elseif example == 2
println("Example 2, sprand - Not working")
  A = sprand(n, n, 0.3) + speye(n,n)
  b = sprand(n, 0.9)

elseif example == 3
println("Example 3, force cast to sparse, full right hand side - Working")
  A = rand(n,n)
  A = sparse(A)

  b = rand(n)
end

println(typeof(A))
println(size(A))
println(typeof(b))
println(size(b))
x = \(A,b)

