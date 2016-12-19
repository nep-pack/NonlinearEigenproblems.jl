function [T, T_prime] = ex1_e_jarl(lambda, additional_flag)
%[T, T_prime] = ex1_e_jarl(lambda, additional_flag)
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
  A0 = [-16, -4, 7
        -14,  7, 13
         6,   8, 7];
  A1 = [0, 0, 0
        0, 0, 0
        0, 0, 0];
  A2 = [ 2, -6,  1
        -2,  22, 11
         7,   -1, 1];
  A3 = [-4,   3,  12
        -17, -11, 0
         1,  -1,  3];

  T = A0 + lambda*A1 + lambda^2 * A2 + lambda^3 * A3;
  T_prime = A1 + 2*lambda*A2 + 3*lambda^2 * A3;
elseif nargin == 2
  T = [3 -2 0.1;
       5  3 0];
  T_prime = -0.4 + 0.6i;
else
  error('Wrong number of inputs.');
end

end
