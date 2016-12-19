function [lambda_k_p_1, v_k_p_1] = inverse_iteration(lambda_k, v_k, d, T, T_prime, ~, ~, ~, ~)
%[lambda_k_p_1, v_k_p_1] = inverse_iteration(lambda_k, v_k, d, T, T_prime, ~, ~)
% Performs one step in the Augmented Newton method for non-linear
% eigenvalue problems. Also known as Inverse Iteration.
%
% lambda_k_p_1  =  is the updated eigenvalue estimate
% v_k_p_1  =  is the updated eigenvector estimate
%
% lambda_k  =  current eigenvalue estimate
% v_k  =  current eigenvector estimate
% d  =  a normalization vector used by the algorithm (d' * v_k = 1)
% T  =  The matrix function T(lambda_k) for which we are searching looking for the eigenvalue
% T_prime  =  The derivative of T with respect to lambda
% ~  =  interface dummy
% ~  =  interface dummy
% ~  =  interface dummy
% ~  =  interface dummy

  first_temp_vec = T_prime * v_k;
  second_temp_vec = T\first_temp_vec;
  
  alpha_inverse = d' * second_temp_vec;
  alpha = 1/alpha_inverse;
  
  v_k_p_1 = second_temp_vec * alpha;
  lambda_k_p_1 = lambda_k - alpha;

end