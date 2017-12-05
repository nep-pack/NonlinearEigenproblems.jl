function y=cheb_polyval(x,c,rho,gamma)
% cheb_polyval: evaluate a polynomial written in 
% shifted-and-scaled Chebyshev basis
% Given the polynomial
% p(x)=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
% this function returns the vector y such that y(i)=p(x(i)) for every 
% component of x

    % Chebyshev polynomials of the first kind
    T=@(n,x) cos(n*acos(x));

    % Chebyshev polynomials (shifted and rescaled) of the first kind
    TT=@(n,x) T(n,rho*x+gamma);
    y=zeros(size(x));
    for i=1:length(x)
        y(i)=0;
        for j=0:length(c)-1
            y(i)=y(i)+c(j+1)*TT(j,x(i));
        end
    end
end