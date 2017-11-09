function [ A0, A1, y0comp ] = generate_problem(  )
%GENERATE_PROBLEM Summary of this function goes here
%   Detailed explanation goes here


n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];
tau=1;



%%%  Compute a very exact solution 
a=-1;
b=0; k=2/(b-a) ; c=(a+b)/(a-b);
cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));

end

