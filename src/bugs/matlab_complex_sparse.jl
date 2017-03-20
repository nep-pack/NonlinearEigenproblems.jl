workspace()

using MATLAB

A = 1im*sprand(10,10,0.15)
B = sprand(10,10,0.15)

aa = mxarray(A)
println(A)
bb = mxarray(B)
println(B)

@mput aa bb
@matlab begin
disp(aa)
disp(bb)

end



C = sprand(10,10,0.15).*(rand()+1im*rand())
println(C)
cc1 = mxarray(real(C))
cc2 = mxarray(imag(C))

@mput cc1 cc2
@matlab begin
cc = cc1 + 1i*cc2
disp(cc)

end
