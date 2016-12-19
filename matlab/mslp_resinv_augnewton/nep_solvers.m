clear;
close all;
clc;

% Set solver and problem
nep_solver_function = str2func('residual_inverse_iteration'); %inverse_iteration residual_inverse_iteration  mslp
function_handle_for_T = str2func('ex2_e_jarl'); %ex1_e_jarl  ex2_e_jarl  ex1_a_neum  ex2_a_neum
verbose = true;


% Run the algorithm on the problem
[A, B] = function_handle_for_T('give the problem number as next input', 2);
lambda = B;
v = (A(1,:))';
d = (A(2,:))';

T_sigma = function_handle_for_T(lambda);
[L, U] = lu(T_sigma);


fprintf('Initial values:\n');
fprintf('lambda = %d + %di\n', real(lambda), imag(lambda));
fprintf('v_k = %d + i*%d\n', real(v(1))', imag(v(1))');
for idx = 2:length(v)
  fprintf('      %d + i*%d\n', real(v(idx)), imag(v(idx)));
end
[T, T_prime] = function_handle_for_T(lambda);
rhs_of_eig_eqn = T*v;
fprintf('|| T(lambfa_k) v_k || = %d\n', norm(rhs_of_eig_eqn));

step_crit = 1e-6;
max_iter_count = 500;
lambda_store = NaN(1,max_iter_count);
conv_crit_1 = false;
conv_crit_2 = false;
v_k_m_1 = inf;
lambda_k_m_1 = inf;
k = 1;
while (k <= max_iter_count) && (~conv_crit_1 || ~conv_crit_2)
  conv_crit_1 = (norm(v-v_k_m_1) < step_crit * norm(v));
  conv_crit_2 = (abs(lambda-lambda_k_m_1) < step_crit * abs(lambda));
  
  d = v/norm(v)^2; %Change the normalization vector
  v_k_m_1 = v;
  lambda_store(k) = lambda;
  lambda_k_m_1 = lambda;
  [T, T_prime] = function_handle_for_T(lambda);
  
  [lambda, v] = nep_solver_function(lambda, v, d, T, T_prime, function_handle_for_T, T_sigma, L, U);
  
  if(verbose)
    fprintf('\nIteration %i of %i \n', k, max_iter_count);
    fprintf('lambda = %d + %di\n', real(lambda), imag(lambda));
    rhs_of_eig_eqn = T*v;
    fprintf('|| T(lambfa_k) v_k || = %d\n', norm(rhs_of_eig_eqn));   
  end
  
  k = k + 1;
end
lambda_store(k) = lambda;

fprintf('v_k = %d + i*%d\n', real(v(1))', imag(v(1))');
for idx = 2:length(v)
  fprintf('      %d + i*%d\n', real(v(idx)), imag(v(idx)));
end


plot(1:k,real(lambda_store(~isnan(lambda_store))), 1:k,imag(lambda_store(~isnan(lambda_store))));
legend('real part', 'imaginary part')