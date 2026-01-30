function [alpha_est_new,rho_est_new,sigma_rho_est_new] = ...
                                                argmaxQ(alpha_est,...
                                                        rho_est,...
                                                        sigma_rho_est, ...
                                                        P, ...
                                                        learning_rate,...
                                                        max_iterations, ...
                                                        tol)
% Calculates the values for alpha_tilde, rho_tilde and sigma_rho_tilde
% that maximize Q(alpha_tilde, rho_tilde, sigma_rho_tilde | alpha_est,
% rho_est, and sigma_rho_est, where Q is defined in the paper. 
%
% alpha_est: a r-dimensional vector of mixing probabilities
% rho_est: a cell array of L cells, where rho_est{l} is a 
%          vector of length k(k-1)/2. (A vectorized correlation matrix)
% sigma_rho_est: a cell array of L cells, where sigma_est{l} is a 
%             (k(k-1)/2) x (k(k-1)/2) covariance matrix.
%
% Projected gradient descent will be used to ensure that rho_tilde remains  
% a vectorization of a valid correlation matrix

% Step 1: Initialize the tilde variables (variables being optimized over)

L = length(rho_est);
r = length(alpha_est);
len_rho = length(rho_est{1}); % = k(k-1)/2

alpha_tilde = alpha_est;  % Initialize to current estimates
rho_tilde = rho_est;
sigma_rho_tilde = sigma_rho_est;

iteration = 1;
norm_grad = Inf;
% Step 2: Loop until max_iterations or ||gradient|| < tol
while iteration < max_iterations & norm_grad >= tol
    iteration = iteration + 1;

    % Step 2.1: Compute derivative of Q for optimization (tilde) variables
                gradient_alpha_tilde = zeros(1,r); % r vector
                gradient_rho_tilde = cell(1,L); % k(k-1)/2-vector for each 
                                                % cell l
                for l = 1:L
                    gradient_rho_tilde{l} = zeros(len_rho,1);
                end
    % Step 2.2: Apply the gradient step to optimization (tilde) variables
                alpha_tilde = alpha_tilde - ...
                             learning_rate*gradient_alpha_tilde;
                for l=1:L
                    rho_tilde{l} = rho_tilde{l} - ...
                                   learning_rate*gradient_rho_tilde{l};
                end
    % Step 2.3: Apply the correlation projection to rho_tilde. Will need
    %           to convert to matrix form, apply projection, and convert
    %           back to vector form. 
                for l=1:L
                    matrix_rho_tilde = vecLInverse(rho_tilde{l});
                    matrix_rho_tilde = ensureValidNearCorrInput(...
                                            matrix_rho_tilde,.01);
                    matrix_rho_tilde = nearcorr(matrix_rho_tilde);
                    rho_tilde{l} = vecL(matrix_rho_tilde);
                end
end
alpha_est_new = alpha_tilde;
rho_est_new = rho_tilde;
sigma_rho_est_new = sigma_rho_tilde;

end
