a=rand()+rand()*im
b=rand()+rand()*im
c=rand()+rand()*im
d=rand()+rand()*im


f=z->(a*z+b)/(c*z+d)
finv=z->(d*z-b)/(-c*z+a)

z=rand()+rand()*im

println(z-f(finv(z)))
println(z-finv(f(z)))
