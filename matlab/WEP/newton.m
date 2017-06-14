function [ x ] = newton( f, fp, x0 )
%NEWTON implements Newton's method
%  f = the function
%  fp = derivative
%  x0 = initial guess


x=x0;
d=1;
tol=1e-14;

i=0;
while (abs(d)>tol) && (i<200)
    
    d=f(x)/(fp(x));
    x=x-d;
    
    i=i+1;
end

end

