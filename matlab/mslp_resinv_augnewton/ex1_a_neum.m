function [T, T_prime] = ex1_a_neum(lambda, additional_flag)
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
% Example comes from A. Neumaier, Residual Inverse Iteraqtion for the
% Nonlinear Eigenvalue Problem, SIAM Journal on NA, Oct. 1985.

if nargin == 1
  A0 = gallery('frank', 11);

  T = A0 - lambda*eye(11);
  T_prime = -eye(11);
elseif nargin == 2
  T = [-1, 0, 1, 0, -1, 0, 1, 0, -1, 0, 1;
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 1];
  T_prime = 0.86;
else
  error('Wrong number of inputs.');
end

end