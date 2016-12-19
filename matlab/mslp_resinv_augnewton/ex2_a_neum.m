function [T, T_prime] = ex2_a_neum(lambda, additional_flag)
%[T, T_prime] = ex2_a_neum(lambda, additional_flag)
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
  lambda2 = lambda^2;
  lower_T = [ 0                      0                      0                    0                      0
              2*lambda2+2*lambda+2   0                      0                    0                      0
             -lambda2+lambda-1       2*lambda2+2*lambda+3   0                    0                      0
              lambda2+2*lambda+2    -2*lambda2+lambda-1    -lambda2-2*lambda+2   0                      0
              3*lambda2+lambda-2    -lambda2+3*lambda-2     lambda2-2*lambda-1   2*lambda2+3*lambda+1   0];
  diagonal_T = diag([-10*lambda2+lambda+10, -11*lambda2+lambda+9, -12*lambda2+10, -10*lambda2+2*lambda+12, -11*lambda2+3*lambda+10]);
  T =  diagonal_T + lower_T + lower_T';

  lower_T_prime = [ 0            0             0            0            0
                    4*lambda+2   0             0            0            0
                   -2*lambda+1   4*lambda+2    0            0            0
                    2*lambda+2   -4*lambda+1  -2*lambda-2   0            0
                    6*lambda+1   -2*lambda+3   2*lambda-2   4*lambda+3   0];
  diagonal_T_prime = diag([-20*lambda+1, -22*lambda+1, -24*lambda, -20*lambda+2, -22*lambda2+3]);
  T_prime =  diagonal_T_prime + lower_T_prime + lower_T_prime';
  
elseif nargin == 2
  switch additional_flag
    case 0
      error('No guess has been found for this value');
    case 1
      T = [-1, 0, 1, 0, -1;
            0, 0, 1, 0,  0];
      T_prime = -1.17;
    case 2
      T = [-1, 0, 1, 1, -1;
            0, 0, 1, 1,  1];
      T_prime = -0.8;
    case 3
      error('No guess has been found for this value');
    case 4
      T = [-1, 0, 1, 0, -1;
            1, 0, 1, 0,  0];
      T_prime = -0.75;
    case 5
      T = [-1, 0, 1, 0, -1;
            1, 0, 0, 0,  0];
      T_prime = 0.4967;
    case 6
      T = [1, 1, 1, 1, 1;
           0, 0, 0, 1, 0];
      T_prime = 0.77;
    case 7
      T = [1, 1, 1, 1, 1;
           0, 0, 0, 0, 1];
      T_prime = 1;
    case 8
      T = [1, -0.5, 1, 0.4, 0.3;
           0.3, 1, 0.2, 0, 0];
      T_prime = 1.6;
    case 9
      T = [1, 1, 1, 1, 1;
           0, 0, 0, 0, 1];
      T_prime = 1.7;
    otherwise
        error('This is an invalid flag.');
  end
else
  error('Wrong number of inputs.');
end

end