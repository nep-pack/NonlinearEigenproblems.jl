function [coeffs, nep_funs, sol]=distributed_delay1(quadv_tol)

    if nargin<1, quadv_tol=eps; end;
    A1=[2.5    2.8   -0.5
        1.8    0.3    0.3
        -2.3   -1.4    3.5];
    A2=[1.7    0.7   -0.3
        -2.4   -2.1   -0.2
        2.0    0.7    0.4];
    A3=[1.4   -1.3    0.4
        1.4    0.7    1.0
        0.6    1.6    1.7];
    
    coeffs={-eye(3),A1,A2,A3};
    nep_funs=@(l) nep_funs_internal(l,quadv_tol);

    if (nargout > 2)
        % load solutions for solution when quadv_tol=eps 
        load distributed_delay1;
        
    end
    
function varargout =nep_funs_internal(lam,quadv_tol)


    lam=lam(:);
    n=length(lam);

    distr_kernel=@(s) exp((s+0.5).^2)-exp(1/4);
    
    distr_vals=quadv(@(s) exp(lam*s)*distr_kernel(s),-1,0,quadv_tol);

    varargout{1}=[lam, ones(n,1), exp(-lam), distr_vals];
        
    if (nargout >= 2),
        distr_vals_der=quadv(@(s) exp(lam*s)*(s*distr_kernel(s)),-1,0,quadv_tol);
        varargout{2}=[ones(n,1),  zeros(n,1), -exp(-lam), distr_vals_der];

    end
    for i=2:nargout-1
        distr_vals_der=quadv(@(s) exp(lam*s)*((s.^i)*distr_kernel(s)),-1,0,quadv_tol);
        varargout{i+1}=[zeros(n,2), ((-1)^i)*exp(-lam), distr_vals_der];
    end

    
    
    