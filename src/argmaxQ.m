function [alpha_est_new,rho_est_new,sigma_rho_est_new] = ...
                                                argmaxQ(alpha_est,...
                                                        rho_est,...
                                                        sigma_rho_est, ...
                                                        X, ...
                                                        P, ...
                                                        w, ...
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

L = length(X);
r = length(alpha_est);
len_rho = length(rho_est{1}); % = k(k-1)/2

% Initialize grad variables to current estimates

rho_tilde = rho_est;
sigma_rho_tilde = sigma_rho_est;

iteration = 1;
norm_grad = Inf;

% Step 2: Calculate the alpha_tilde that maximizes the alpha_tilde terms

% Calculate Prob( X | Gamma_lj = 1, theta). This is a Lxr matrix
% This is normal with mean mu = (P_l rho_j)_{1...len(X)} and Sigma = 
% (P_l Sigma_rho,j P_l')_{1...len(X),1...len(X)}, that is the XX block
% of the joint (X,Z) distribution 
pr_x_given_g = zeros(L,r);

% THIS LOOP IS SLOW: FIX 
mu_z_given_x = cell(L,r);
for j=1:r
    for l=1:L
        x_len = length(X{l});
        mu_est = P{l}*rho_est{j};
        mu_X_est = mu_est(1:x_len); 
        mu_Z_est = mu_est((x_len+1):end);
        sig_est = w(l)*P{l}*sigma_rho_est{j}*P{l}'; 
        sig_XX_est = sig_est(1:x_len,1:x_len);
        sig_ZX_est = sig_est(1:x_len,(x_len+1):end);
        % Calculate conditional mean Z|X, which we will need for rho_tilde
        mu_Z_given_X{l,j} = ...
                mu_Z_est + sig_ZX_est'*...
                ((sig_XX_est)\(X{l}-mu_X_est)); %mu_{Z given X} or mu_{2|1}
        % Calculate likelihood
        pr_x_given_g(l,:) = mvnpdf(X{l},mu_X_est,sig_XX_est)+eps;
    end
end

% pr_g_given_x: Prob (Gamma_lj = 1 | X, theta). Calculate using Bayes
% theorem
% pr_g_given_x_lsum: Sum_{l=1}^L Prob (Gamma_lj = 1 | X, theta)
% pr_g_given_x_ljsum Sum_{j=1}^r Sum_{l=1}^L Prob (Gamma_lj = 1 | X, theta)
% In the following, pr_x_given_g is Lxr and alpha_est is an r-vector
% Each column of pr_x_given_g is multiplied by corresponding scalar element
% of alpha_est (jth) and then each row is divided by corresponding scalar
% element of the L-vector sum(pr_x_given_g.*alpha_est,2); 
pr_g_given_x = (pr_x_given_g.*alpha_est)./sum(pr_x_given_g.*alpha_est,2); 
pr_g_given_x_lsum = sum(pr_g_given_x,1); % sum over l index
pr_g_given_x_ljsum = L; % sum over j=1,..,r is just 1, so l,j sum is L

% See expression for multinomial MLE 
alpha_tilde = pr_g_given_x_lsum/pr_g_given_x_ljsum;

% Step 3: Loop until max_iterations or ||gradient|| < tol
while iteration <= max_iterations & norm_grad >= tol
    iteration = iteration + 1;

    % Step 3.1: Compute derivative of Q for remaining optimization (tilde) 
    %           variables (rho_tilde)
    
                gradient_rho_tilde = cell(1,r); % k(k-1)/2-vector for each 
                                                % cell l
                for i = 1:r
                    % Compute the sum of the l-summands
                    total = 0;
                    for l = 1:L
                        % l'th term in rho gradient sum 
                        total = total + ...
                           pr_g_given_x(l,j)*...
                           P{l}'*inv(w(l)*P{l}*sigma_rho_tilde{i}*P{l}')...
                           *([X{l} ; mu_Z_given_X{l,i}] - ...
                           P{l}*rho_tilde{i});
                    end
                    gradient_rho_tilde{i} = -total;%zeros(len_rho,1);
                end

    % Step 3.2: Apply the gradient step to remaining optimization (tilde) 
    %           variables (rho_tilde)

                for j=1:r
                    rho_tilde{j} = rho_tilde{j} - ...
                                   learning_rate*gradient_rho_tilde{j};
                end

    % Step 3.3: Apply the correlation projection to rho_tilde. Will need
    %           to convert to matrix form, apply projection, and convert
    %           back to vector form. 

                for j=1:r
                    matrix_rho_tilde = vecLInverse(rho_tilde{j});
                    matrix_rho_tilde = ensureValidNearCorrInput(...
                                            matrix_rho_tilde,.01);
                    matrix_rho_tilde = nearcorr(matrix_rho_tilde);
                    rho_tilde{j} = vecL(matrix_rho_tilde);
                end

end

% Step 4: Return the new estimates
alpha_est_new = alpha_tilde;
rho_est_new = rho_tilde;
sigma_rho_est_new = sigma_rho_tilde;

end
