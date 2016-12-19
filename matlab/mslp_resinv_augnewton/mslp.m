function [lambda_k_p_1, v_k_p_1] = mslp(lambda_k, ~, ~, T, T_prime, ~, ~, ~, ~)
%[lambda_k_p_1, v_k_p_1] = mslp(lambda_k, ~, ~, T, T_prime, ~, ~)
% Performs one step in the Augmented Newton method for non-linear
% eigenvalue problems. Also known as Inverse Iteration.
%
% lambda_k_p_1  =  is the updated eigenvalue estimate
% v_k_p_1  =  is the updated eigenvector estimate
%
% lambda_k  =  current eigenvalue estimate
% ~  =  interface dummy
% ~  =  interface dummy
% T  =  The matrix function T(lambda_k) for which we are searching looking for the eigenvalue
% T_prime  =  The derivative of T with respect to lambda
% ~  =  interface dummy
% ~  =  interface dummy
% ~  =  interface dummy
% ~  =  interface dummy

%   normal_linear_eig_matrix = T_prime \ T;
%   [V, eig_vals] = eig(normal_linear_eig_matrix, 'vector');
  [V, eig_vals] = eig(T, T_prime, 'vector');
  [~, idx] = min(abs(eig_vals));
  
  v_k_p_1 = V(:,idx);
  mu = eig_vals(idx);
  lambda_k_p_1 = lambda_k - mu;
end