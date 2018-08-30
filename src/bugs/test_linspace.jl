using MATLAB

if !isdefined(:nx)
    nx = 231
end
if !isdefined(:delta)
    delta = 0.1
end

@mput nx delta
@matlab begin
    xp = 2/pi+0.4;
    xp = xp + delta;
    xm = 0;
    xm = xm - delta;
    X_m = linspace(xm, xp, nx+2)';
    hhh_m = X_m(3)-X_m(2);
    ccc_m = 1/hhh_m^2
@matlab end
@mget hhh_m ccc_m X_m

temp = 2/pi +0.4 + 0.1;
X = range(-0.1, stop = temp, length = nx+2)
hhh = step(X)
X = collect(X)
hhh_naive = X[3] - X[2];
ccc = 1/hhh^2

println("Matlab step size: ", hhh_m)
println("Julia step size: ", hhh)
println("Julia (naive) step size: ", hhh_naive)
println("Difference: ", hhh-hhh_m)
println("Relative difference: ", abs(hhh-hhh_m)/hhh)
println("Difference (naive): ", hhh_naive-hhh_m)
println("Relative difference (naive): ", abs(hhh_naive-hhh_m)/hhh_naive)
println("")
println("Matlab 2nd derivative: ", ccc_m)
println("Julia 2nd derivative: ", ccc)
println("Difference: ", ccc-ccc_m)
println("Relative difference: ", abs(ccc-ccc_m)/ccc)
println("")
println("The difference in points:")
println(X-X_m)
