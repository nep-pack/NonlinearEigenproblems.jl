# This file is a direct conversion of the work of
# Steinlechner described in
#    Steinlechner, A boundary element method for solving PDE eigenvalue problems, M. Steinlechner, bachelor thesis, ETH Zürich, 2010, http://sma.epfl.ch/~anchpcommon/students/steinlechner.pdf
# and used in
# "Chebyshev interpolation for nonlinear eigenvalue problems",
# Effenberger and Kressner
# BIT Numerical Mathematics, December 2012, Volume 52, Issue 4, pp 933–951
#
# https://anchp.epfl.ch/wp-content/uploads/2018/05/chebapprox.tar.gz
#

function solidAngle(P1,P2,P3)
    numer = abs.(P1[1, :].*P2[2, :].*P3[3, :] - P1[1, :].*P2[3, :].*P3[2, :]+
                 P1[2, :].*P2[3, :].*P3[1, :] - P1[2, :].*P2[1, :].*P3[3, :]+
                 P1[3, :].*P2[1, :].*P3[2, :] - P1[3, :].*P2[2, :].*P3[1, :]);

    lP1 = sqrt.(sum(P1.^2, dims=1));
    lP2 = sqrt.(sum(P2.^2, dims=1));
    lP3 = sqrt.(sum(P3.^2, dims=1));

    denom = lP1.*lP2.*lP3 + lP1.*sum(P2.*P3, dims=1) +
        lP2.*sum(P1.*P3, dims=1) + lP3.*sum(P1.*P2, dims=1);
    solang = 2*atan.(vec(numer), vec(denom));
    idx = (solang .< 0); # Fix branch
    solang[idx] = solang[idx] .+ 2*pi;
    return solang;
end


function  deHoop(x, tri)

    m = size(x, 2);
    R1 = tri.P1[:]*ones(Float64,1, m) - x;
    R2 = tri.P2[:]*ones(Float64,1, m) - x;
    R3 = tri.P3[:]*ones(Float64,1, m) - x;

    normR1 = sqrt.(sum(R1.^2, dims=1));
    normR2 = sqrt.(sum(R2.^2, dims=1));
    normR3 = sqrt.(sum(R3.^2, dims=1));

    dist = vec(abs.(tri.aux.normal' * R1));

    solang = solidAngle(R1, R2, R3);

    dot_R1_Nu1 =  tri.aux.nu1' * R1;
    dot_R2_Nu2 =  tri.aux.nu2' * R2;
    dot_R3_Nu3 =  tri.aux.nu3' * R3;

    dot_R1_Tau1 = tri.aux.tau1' * R1;
    dot_R2_Tau2 = tri.aux.tau2' * R2;
    dot_R3_Tau3 = tri.aux.tau3' * R3;

    dot_R2_Tau1 = tri.aux.tau1' * R2;
    dot_R3_Tau2 = tri.aux.tau2' * R3;
    dot_R1_Tau3 = tri.aux.tau3' * R1;

    # The deHoop formula
    F0=-dist .* solang;
    F1=vec(dot_R1_Nu1 .* log.( (normR2 + dot_R2_Tau1) ./ (normR1 + dot_R1_Tau1) ));
    F2=vec(dot_R2_Nu2 .* log.( (normR3 + dot_R3_Tau2) ./ (normR2 + dot_R2_Tau2) ));
    F3=vec(dot_R3_Nu3 .* log.( (normR1 + dot_R1_Tau3) ./ (normR3 + dot_R3_Tau3) ));
    F = F0+F1+F2+F3

    return F
end
