function [ err ] = backward_error( nep, lambda, v )
%BACKWARD_ERROR computes an error estimate for the NEP


    sM=nep.sM(lambda)-nep.d0; sM=abs(sM); sM=sum(sM);
    sP=nep.sP(lambda)-nep.d0; sP=abs(sP); sP=sum(sP);
    
    norm_1_A0 = 4*1/nep.hz^2 + 4*1/nep.hx^2 + max(max(abs(nep.K)));
    norm_1_A1 = 2*1/nep.hz;
    norm_1_A2 = 1;
    
    alpha=norm_1_A0 + abs(lambda)*norm_1_A1 + abs(lambda^2)*norm_1_A2 +...
          norm(nep.C1,1) + norm(nep.C2T,1) + sM + sP + 2*abs(nep.d0);

    err=norm(nep.M(lambda,v))/(alpha*norm(v));


end

