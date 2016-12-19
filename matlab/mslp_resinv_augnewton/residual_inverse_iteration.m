function [lambda_k_p_1, v_k_p_1] = residual_inverse_iteration(lambda_k, v_k, d, ~, ~, function_handle_for_T, T_sigma, L, U)
%residual_inverse_iteration(lambda_k, v_k, d, ~, ~, function_handle_for_T, A_sigma)
% Performs one step in the Residual Inverse Iteration.
%
% lambda_k_p_1  =  is the updated eigenvalue estimate
% v_k_p_1  =  is the updated eigenvector estimate
%
% lambda_k  =  current eigenvalue estimate
% v_k  =  current eigenvector estimate
% d  =  a normalization vector used by the algorithm (d' * v_k = 1)
% ~  =  interface dummy
% ~  =  interface dummy
% function_handle_for_T  =  Function handle to be able to genereate T and
%                           it's derivative T_prime
% T_sigma  =  An matrix T, but with a fixed shift sigma instead of updating
%             it with lambda.
% L  =  L form an LU-factorization of a fixed-shift T
% U  =  U form an LU-factorization of a fixed-shift T

  lambda_k_p_1 = lambda_k;
  h = lambda_k;
  max_nrof_iters = 50;
  n = 1;
  %Solve: d' * A_sigma^(-1) * A * v_k = 0
  while norm(h)>1.e-6*norm(lambda_k_p_1) && n < max_nrof_iters
    [T, T_prime] = function_handle_for_T(lambda_k_p_1);
    left_vector = d' / T_sigma;
    %function = 0
    f = left_vector * T * v_k;
    %Jacobian
    J = left_vector * T_prime * v_k;

    h = -J \ f;
    lambda_k_p_1 = lambda_k_p_1 + h;

    n = n+1;
  end
  if( norm(h,inf)>1.e-6*norm(lambda_k_p_1,inf) && n >= max_nrof_iters)
    error('Did not converge');
  end

  [T, ~] = function_handle_for_T(lambda_k_p_1);
  r = T*v_k;
  dv = U\(L \ r);
  v_k_p_1 = v_k - dv;
  v_k_p_1 = v_k_p_1 / (d' * v_k_p_1);
  
end