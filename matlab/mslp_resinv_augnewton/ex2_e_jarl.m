function [T, T_prime] = ex2_e_jarl(lambda, additional_flag)
%[T, T_prime] = ex1_a_neum(lambda, additional_flag)
%
% If only lambda is given, this function calculates the T-matrix and it's
% derivative, T and T_prime respecitvely.
%
% If two arguments are given the function will retun initial guesses to
% start an NEP-solver.
% T will contain two rows, first one is an initial guess to the
% eigenvector, second one is a normailzation vector d used in some methods.
% T_prime  will contain an initial guess for the eigenvalue.
%
% Example comes from E. Jarlebring, Convergence factors of Netwon methods
% for nonlinear eigenvalue problems, Linear Algebra and its Applications,
% 2012

if nargin == 1
  denom = 8+5*pi;
  a1 = 2/5 *(65*pi + 32)/(denom);
  a2 = 9*pi^2*(13+5*pi)/(denom);
  a3 = 324/5 *pi^2*(5*pi+4)/(denom);

  b1 = (260*pi + 128 + 225*pi^2)/(10*denom);
  b2 = 45*pi^2/denom;
  b3 = 81*pi^2*(40*pi + 32 + 25*pi^2)/(10*denom);

  A0 = [ 0    1    0
         0    0    1
        -a3, -a2, -a1];
  A1 = [ 0,   0,   0
         0,   0,   0
        -b3, -b2, -b1];


  T = -lambda*eye(3) + A0 + A1*exp(-lambda);
  T_prime = -eye(3) - A1*exp(-lambda);

elseif nargin == 2
  
  T = [5+1i -3 05+3i;
       1 1 1];
  T_prime = 0.4 + 8.7i;
else
  error('Wrong number of inputs.');
end

end
