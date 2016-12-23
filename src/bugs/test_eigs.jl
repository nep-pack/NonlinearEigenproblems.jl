# TEST: eigs with JULIA and with MATLAB
# with JULIA eigs does not work
# see also https://github.com/JuliaInterop/MATLAB.jl

workspace()

a = sprand(100,100,0.3);
b = sprand(100,100,0.3);
d,v = eigs(a,b);

using MATLAB
aa=mxarray(a)
bb=mxarray(b)
@mput aa bb;
@matlab (vv,dd)=eigs(aa,bb); 
@mget dd vv;

println("residual=",norm(a*v[:,2] - d[2]*b*v[:,2]))
println("norm eigenvectors=",norm(v))

println("residual matlab=",norm(a*vv[:,2] - dd[2,2]*b*vv[:,2]))
println("norm eigenvectors matlab=",norm(vv))
