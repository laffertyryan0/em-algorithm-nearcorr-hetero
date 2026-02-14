%% Set path
addpath('./src')

%% Simulate data

tic

% We have k metabolites
k = 20;
% We have L labs 
L = 4;

average_fraction_missing_metabolites = .1; % approx what proportion of 
                                           % metabolites
                                           % are missing (1 = all missing
                                           % 0 = none missing)
                                           % sampled as iid bernoulli

% Each lab has a state vector Gamma_l which is a one hot binary vector
% Gamma_l,j = 1 iff lab is in state j. Let's say j = 1...r, maybe r=2 
r = 2;
alpha = ones(1,r)/r; % The alpha_j are the probabilities of a lab being 
                   % selected as state j


Gamma = mnrnd(1,alpha,L); % multinomial dist with n=1 is categorical dist
                          % Gamma is Lxr with each row a One Hot Encoded
                          % lab state


% Number of subjects in lab l is n_subj(l)
n_subjects = ones(L,1)*100;  % In this case assume 100 per lab

% Each lab recruits n_subjects(l) patients, all with the same state gamma_l
% Why? because gamma_l is specific to a (say) regional population,
% and suppose any two individuals who could be recruited by the same
% lab share this regional characteristic gamma_l (e.g. regional
% diet/climate)

% Imagine that if this study were replicated, the states would be shuffled
% again randomly. Alternatively, one could imagine a fixed unknown state 
% assumption

% Next, for each lab, simulate the subject level data
% Fix a correlation matrix for each state
% To keep things simple, consider all the data to be mean 0, var 1
% It shouldn't affect the results since we are only looking at correlations

rho_state = {}; % cell array where rho_state{j} is the kxk rho matrix for 
                % the jth state

for j=1:r
    rho_state{j} = randomCorrelationMatrix(k); % k = number of metabolites
    assert(min(eig(rho_state{j}))>= 0, ...
        "Non-PSD matrix found for simulation rho_i.")
end
clear j

% Generate subject level data for each lab
subject_data = {}; % subject_data{i} is design matrix for 
                         % ith lab
% Define a mask for each lab: 1 for non-missing, 0 for missing
% subject_data_mask{l} is a binary vector of length n_subjects(l)
% Lab l only "really" collects data for unmasked (mask=1) metabolites
subject_data_mask = {};  % 1 = non-missing, 0 = missing

% Also for each lab, calculate the spearman rank correlation estimate
% The spearman correlation estimates only are reported to the analyst
% We also need a mask for spearman correlation entries, corresponding
% to pairs 1<=(j1,j2)<=k where mask=1 iff metabolites j1,j2 both 
% collected at lab l
% reported_spearman{l}(j1,j2) at a mask=0 location will be 0 and is 
% not meaningful

reported_spearman = {};
reported_spearman_mask = {}; % 1 = non-missing, 0 = missing

for l=1:L
    mu = zeros(1,k);                   % Centered mean for simplicity
    lab_state = find(Gamma(l,:));      % Latent state for this lab = 
                                       % index of the 1 column for l'th lab
    Sigma = rho_state{lab_state};      % Variances all one, so Sigma = rho
    subject_data{l} = mvnrnd(mu,...
                             Sigma, ...
                             n_subjects(l));


    % Generate a random mask matrix for lab l
    p = 1-average_fraction_missing_metabolites; %p = 0 means none are missing
    covariate_mask = rand(1,k) < p; % 1 if non-missing
    subject_data_mask{l} = covariate_mask;
    
    % Mask out the missing data: that is, set them to 0. These 0's
    % are meaningless placeholders (don't interpret as data = 0!)
    subject_data{l}(:,covariate_mask==0)=0; 
    
    % The following has some entries masked out. Those entries are zero
    % and are meaningless
    reported_spearman{l} = corr(subject_data{l},'Type','Spearman');
    reported_spearman_mask{l} = ones(k,k);
    reported_spearman_mask{l}(:,covariate_mask==0)=0;
    reported_spearman_mask{l}(covariate_mask==0,:)=0;
    reported_spearman{l}(reported_spearman_mask{l}==0)=0; % zero out 
                                                          % masked entries

end
clear mu lab_state Sigma l mask p covariate_mask;

% The observable data (available to the analyst) consists of
% reported_spearman, and reported_spearman_mask

data_generation_time = toc;
fprintf("Simulated data generation complete. Took %.2f seconds.\n",...
    data_generation_time);


%% Using simulated data (lab aggregates only), estimate correlation matrix

% We have access to lab sample sizes, reported_spearman and 
% reported_spearman mask in this section. 

% First we need to place the matrices into vector form 
% The vech function flattens the upper triangular entries of a matrix
% and is a bijection from symmetric matrices to vectors

reported_spearman_vecL = cellfun(@vecL,reported_spearman,...
                                    'UniformOutput',false);
reported_spearman_mask_vecL = cellfun(@vecL,reported_spearman_mask,...
                                        'UniformOutput',false);

% For each l, find the permutation matrix P_l that puts missing entries 
% below observed entries. Let [X Z]' = P_l * vecL(reported_spearman)
P = {}; 
X = {}; % Consists of all observed entries. 
        % Z here would be represented as all 0's due to the masking, but 
        % really Z means unobserved data. 

for l=1:L
    P{l} = getMaskOrderingMatrix(reported_spearman_mask_vecL{l});
    num_observed = sum(reported_spearman_mask_vecL{l});
    X_Z = P{l}*reported_spearman_vecL{l};
    X{l} = X_Z(1:num_observed);
end

MAX_EM_ITERATIONS = 300; % Outer loop
MAX_GD_ITERATIONS = 100; % Inner PGD loop
GD_TOLERANCE = .001;
GD_LEARNING_RATE = 0.1;

% Initialize EM parameters
alpha_est = rand(1,r);
alpha_est = alpha_est/sum(alpha_est); % Random initialization
rho_est = cell(1,r);
sigma_rho_est = cell(1,r);
for j = 1:r
    rho_est{j} = vecL(randomCorrelationMatrix(k)); % Random initialization
    sigma_rho_est{j} = speye(k*(k-1)/2);
end

w = sqrt(mean(n_subjects))./sqrt(n_subjects); % Lab-wise weighting factor for variances (L vector)


for em_iter=1:MAX_EM_ITERATIONS
    % Update the EM using the following formula: 
    % theta^{(t+1)} = argmax_{theta_tilde} Q(theta_tilde | theta^{(t)})
    % Here, theta represents the current estimate of parameters:
    % alpha_est: a r-dimensional vector of mixing probabilities
    % rho_est: a cell array of r cells, where rho_est{j} is a 
    %          vector of length k(k-1)/2. (A vectorized correlation matrix)
    % sigma_rho_est: a cell array of r cells, where sigma_est{j} is a 
    %              (k(k-1)/2) x (k(k-1)/2) covariance matrix.
    % For now, we will not estimate sigma_rho_est, but just let it remain 
    % fixed. Later, we will try to show that this is a good approximation.
    % X: The observed data (a cell array where each lab is a cell)
    % P: The cell array of permutation matrices that rearrange X and Z
    % w: lab-wise weighting that multiplies the covariance matrices
    % GD_LEARNING_RATE: The rate parameter for gradient descent
    % MAX_GD_ITERATIONS: Number of iterations for GD in each EM step. This
    %                    can be small, since it is not mandatory for the
    %                    gradient descent to converge in each step
    % GD_TOLERANCE: If the gradient steps fall below this tolerance, the 
    %               gradient descent loop will cease.

    [alpha_est,rho_est,sigma_rho_est] = argmaxQ(alpha_est,...
                                                rho_est,...
                                                sigma_rho_est, ...
                                                X, ...
                                                P, ...
                                                w, ...
                                                GD_LEARNING_RATE,...
                                                MAX_GD_ITERATIONS, ...
                                                GD_TOLERANCE);
    fprintf("EM Iteration Number: %d\n",em_iter);
    fprintf("Current alpha estimate: ");
    disp(alpha_est);
    fprintf("\n");
end
