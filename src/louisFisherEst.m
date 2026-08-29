function rho_FI_est = ...
                                                louisFisherEst(alpha_est,...
                                                        rho_est,...
                                                        sigma_rho_est, ...
                                                        X, ...
                                                        P, ...
                                                        w)
% Estimate the observed Fisher information using the Louis trick
% Should be used only after EM convergence
% alpha_est: a r-dimensional vector of mixing probabilities
% rho_est: a cell array of L cells, where rho_est{l} is a 
%          vector of length k(k-1)/2. (A vectorized correlation matrix)
% sigma_rho_est: a cell array of L cells, where sigma_est{l} is a 
%             (k(k-1)/2) x (k(k-1)/2) covariance matrix.
%
% NO GRAD
% Step 1: Initialize the tilde variables (variables being optimized over)

L = length(X);
r = length(alpha_est);
len_rho = length(rho_est{1}); % = k(k-1)/2
k = floor((1+sqrt(1+8*len_rho))/2);


% DEBUG %%%%
 % if mod(em_iteration,100)==0
 %    for j = 1:r
 %        sigma_rho_tilde{j} = (1/mean(w))*getCovarianceOfCorrelations(...
 %                                vecLInverse(rho_est{j}),1/mean(w));
 %    end
 % end
%%%%%%%%%%%%


iteration = 1;
norm_grad = Inf;

% Step 2: Calculate the alpha_tilde that maximizes the alpha_tilde terms

% Calculate Prob( X | Gamma_lj = 1, theta). This is a Lxr matrix
% This is normal with mean mu = (P_l rho_j)_{1...len(X)} and Sigma = 
% (P_l Sigma_rho,j P_l')_{1...len(X),1...len(X)}, that is the XX block
% of the joint (X,Z) distribution 
pr_x_given_g_abs = zeros(L,r);
pr_x_given_g_rel = zeros(L,r);
pr_x_given_g_log = zeros(L,r);


% THIS LOOP IS SLOW: FIX 
mu_Z_given_X = cell(L,r);
sig_Z_given_X = cell(L,r);
sig_XZ_given_X = cell(L,r);
for l=1:L
    % The number of non-missing entries
    x_len = length(X{l});
    % Calculate the mu and sig (permuted vectors) with appropriate
    % partitions for the first mixing component. To be used for 
    % denominators of relative probability ratios
    mu_est_j1 = P{l}*rho_est{1};
    mu_X_est_j1 = mu_est_j1(1:x_len); % Just the non-missing part
    sig_est_j1 = w(l)*P{l}*sigma_rho_est{1}*P{l}'; 
    sig_XX_est_j1 = sig_est_j1(1:x_len,1:x_len); % Just the upper left block
    for j=1:r
        mu_est = P{l}*rho_est{j};
        mu_X_est = mu_est(1:x_len); % Just the non-missing part
        mu_Z_est = mu_est((x_len+1):end); % Just the missing part
        sig_est = w(l)*P{l}*sigma_rho_est{j}*P{l}'; 
        sig_XX_est = sig_est(1:x_len,1:x_len); % Just the upper left block
        sig_XZ_est = sig_est(1:x_len,(x_len+1):end); % Upper right block
        sig_ZZ_est = sig_est((x_len+1):end,(x_len+1):end);
        % Calculate conditional mean Z|X, which we will need for rho_tilde
        mu_Z_given_X{l,j} = ...
                mu_Z_est + sig_XZ_est'*...
                ((sig_XX_est)\(X{l}-mu_X_est)); %mu_{Z given X} or mu_{2|1}
        % Conditional variance
        sig_Z_given_X{l,j} = ...
            (sig_ZZ_est - sig_XZ_est'*inv(sig_XX_est)*...
                                            sig_XZ_est);
        sig_XZ_given_X{l,j} = ...
            [zeros(size(sig_XX_est)) zeros(size(sig_XZ_est));
             zeros(size(sig_XZ_est')) sig_Z_given_X{l,j}];
        % Calculate relative pdf
        pr_x_given_g_rel(l,j) = det(sig_XX_est_j1\sig_XX_est)^(-1/2)*...
                                exp(-(1/2)*(X{l}-mu_X_est)'*...
                                ((sig_XX_est)\...
                                (X{l}-mu_X_est))+...
                                (1/2)*(X{l}-mu_X_est_j1)'*...
                                ((sig_XX_est_j1)\...
                                (X{l}-mu_X_est_j1)));
        %Calculate absolute pdf
        pr_x_given_g_abs(l,j) = (2*pi)^(-x_len/2)*det(sig_XX_est)^(-1/2)*...
                                exp(-(1/2)*(X{l}-mu_X_est)'*...
                                ((sig_XX_est)\...
                                (X{l}-mu_X_est)));
        %Calculate log pdf
        % pr_x_given_g_log(l,j) = log((2*pi)^(-x_len/2)*det(sig_XX_est)^(-1/2))-...
        %                         (1/2)*(X{l}-mu_X_est)'*...
        %                         ((sig_XX_est)\...
        %                         (X{l}-mu_X_est));
    end
end

pr_x_given_g_rel(isinf(pr_x_given_g_rel)) = 1e20;

% pr_g_given_x: Prob (Gamma_lj = 1 | X, theta). Calculate using Bayes
% theorem
% pr_g_given_x_lsum: Sum_{l=1}^L Prob (Gamma_lj = 1 | X, theta)
% pr_g_given_x_ljsum Sum_{j=1}^r Sum_{l=1}^L Prob (Gamma_lj = 1 | X, theta)
% In the following, pr_x_given_g is Lxr and alpha_est is an r-vector
% Each column of pr_x_given_g is multiplied by corresponding scalar element
% of alpha_est (jth) and then each row is divided by corresponding scalar
% element of the L-vector sum(pr_x_given_g.*alpha_est,2); 
% if ~sum(pr_x_given_g_rel>1e20,'all')>0
    pr_g_given_x = (pr_x_given_g_rel.*alpha_est)./...
                    sum(pr_x_given_g_rel.*alpha_est,2);
% else
%     pr_g_given_x = (pr_x_given_g_abs.*alpha_est)./...
%                     sum(pr_x_given_g_abs.*alpha_est,2); 
% end



pr_g_given_x_lsum = sum(pr_g_given_x,1); % sum over l index
pr_g_given_x_ljsum = L; % sum over j=1,..,r is just 1, so l,j sum is L

% See expression for multinomial MLE 
alpha_tilde = pr_g_given_x_lsum/pr_g_given_x_ljsum;

% for debugging, lets compute the actual estimates for a simple mixture
% for j=1:r
%     total = 0;
%     for l=1:L
%         total = total + pr_g_given_x(l,j)*X{l};
%     end
%     %rho_tilde{j} = total./pr_g_given_x_lsum(j);
%     rho_tilde_calculated{j} = total./pr_g_given_x_lsum(j);
% 
% end

% for debugging, lets compute the VARIANCE estimates for a simple mixture

avg_num_samples = 1/mean(w);
% for j=1:r
%     total = 0;
%     for l=1:L
%         total = total + pr_g_given_x(l,j)*(X{l}-rho_tilde{j})^2;
%     end
%     sigma_rho_tilde{j} = avg_num_samples*total./pr_g_given_x_lsum(j);
% end

if 1
end

% rhos1 = [];
% rhos2 = [];

% 
% % Step 3: Loop until max_iterations or ||gradient|| < tol

while iteration <= 0 
    iteration = iteration + 1;

    % Step 3.1: Minimize Q with respect to rho_tilde, not respecting the
    %           constraint
            
                for j=1:r
                    num = 0;
                    denom = 0;
                    for l=1:L
                        num = num + ...
                            pr_g_given_x(l,j)*(1/w(l))*P{l}'*...
                            [X{l} ; mu_Z_given_X{l,j}];
                        denom = denom + ...
                            pr_g_given_x(l,j)*(1/w(l));
                    end
                    rho_tilde{j} = (num)/(denom+eps);
                end


    % Step 3.2: Apply the correlation projection to rho_tilde. Will need
    %           to convert to matrix form, apply projection, and convert
    %           back to vector form. 

                if(NEARCORR_PROJ)
                    for j=1:r
                        matrix_rho_tilde = vecLInverse(rho_tilde{j});
                        % Inverse fisher transformation
                        matrix_rho_tilde = tanh(matrix_rho_tilde); 
                        matrix_rho_tilde(logical(eye(k))) = 1;
                        matrix_rho_tilde = ensureValidNearCorrInput(...
                                                matrix_rho_tilde,.01);
                        matrix_rho_tilde = nearcorr(matrix_rho_tilde);
                        % Fisher transformation
                        matrix_rho_tilde = atanh(matrix_rho_tilde);
                        rho_tilde{j} = vecL(matrix_rho_tilde);
                    end
                end


end


% shouldnt there actually be r different fishers for rho? whats going on
% checked --> its not so bad. look at updated notes. you mostly just
% take out the sums over j

EDel2alphL = zeros(r-1,r-1);
EDel2rhoL = zeros(len_rho,len_rho); % what size? not r? its going to be messed up because stuff is flattened
% so check below when you can. but actually it might be fine since rho is already flattened in the 
% likelihood expression. just need d = r*(r-1)/2? 
EDelalphDelalphL = zeros(r-1,r-1);
EDelrhoDelrhoL = zeros(len_rho,len_rho);
EDel2alphrhoL = zeros(r-1,len_rho); % TODO
EDelalphDelrhoL = zeros(r-1,len_rho); % TODO

for l=1:L
    for i1 = 1:(r-1)
        delta = eye(r-1,r-1);
        EDel2alphL = EDel2alphL - delta(i1,:)*(1/alpha_est(i1)^2) * pr_g_given_x(l,i1);
    end
    for j = 1:r
        EUlj = -pr_g_given_x(l,j) *P{l}'*inv(w(l)*P{l}*sigma_rho_est{j}*P{l}');
        EUlj_given_G1 = ...
            -P{l}'*inv(w(l)*P{l}*sigma_rho_est{j}*P{l}');
        EDel2rhoL = EDel2rhoL -  EUlj*P{l}; 
        EDelrhoDelrhoL = EDelrhoDelrhoL + pr_g_given_x(l,j) * ...
                EUlj_given_G1* ...
                sig_XZ_given_X{l,j}*...
                EUlj_given_G1';
                        
    end
end

expected_hessian = [EDel2alphL    EDel2alphrhoL;
                    EDel2alphrhoL' EDel2rhoL];
expected_squaregrad = [EDelalphDelalphL    EDelalphDelrhoL;
                       EDelalphDelrhoL'    EDelrhoDelrhoL];

% % Step 4: Return the new estimates
% alpha_est_new = alpha_tilde;
% rho_est_new = rho_tilde;
% sigma_rho_est_new = sigma_rho_tilde;

rho_FI_est = expected_hessian - expected_squaregrad;


end
