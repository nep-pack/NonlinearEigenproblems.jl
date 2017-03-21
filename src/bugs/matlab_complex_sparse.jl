workspace()

using MATLAB

A = 1im*sprand(10,10,0.15)
B = sprand(10,10,0.15)

println("One Real and one Complex matrix")

aa = mxarray(A)
println(A)
bb = mxarray(B)
println(B)

@mput aa bb
@matlab begin
disp(aa)
disp(bb)

end

println("Potential work-around")

C = sprand(10,10,0.15).*(rand()+1im*rand())
println(C)
cc1 = mxarray(real(C))
cc2 = mxarray(imag(C))

@mput cc1 cc2
@matlab begin
cc = cc1 + 1i*cc2
disp(cc)

end

println("Works for dense matrices")

D = rand(10,10)+1im*rand(10,10)
println(D)
dd = mxarray(D)

@mput dd
@matlab begin
disp(dd)

end

println("Does not work for BigFloats")

E = [BigFloat(1) 2; 3 4]
println(E)
ee = mxarray(E)

@mput ee
@matlab begin
disp(ee)

end
