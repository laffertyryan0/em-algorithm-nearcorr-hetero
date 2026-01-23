%% Simulate data

% We have k metabolites
k = 200;
% We have L labs 
L = 300;

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

% Generate subject level data for each lab
subject_data = {}; % subject_data{i} is design matrix for 
                         % ith lab
for l=1:L
    mu = zeros(1,k);                   % Centered mean for simplicity
    lab_state = find(Gamma(l,:));      % Latent state for this lab = 
                                       % index of the 1 column for l'th lab
    Sigma = rho_state{lab_state};      % Variances all one, so Sigma = rho
    subject_data{l} = mvnrnd(mu,...
                             Sigma, ...
                             n_subjects(l));
end

%% Using simulated data (lab aggregates only), estimate correlation matrix