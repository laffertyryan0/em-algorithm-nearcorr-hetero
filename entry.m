%% Simulate data

% We have k metabolites
k = 200;
% We have L labs 
L = 300;
% Each lab has a state vector Gamma_l which is a one hot binary vector
% Gamma_l,j = 1 iff lab is in state j. Let's say j = 1...r, maybe r=2 
r = 2;
Gamma = zeros(L,r);
alpha = ones(r,1)/r; % probability of each state
randmemb = randsample(1:r,L,true,alpha)'; % Choose a random member 
                       % of each row
                       % That is for each lab choose a state

idx = sub2ind(size(Gamma),(1:L)',randmemb); % Create indices for assignment
Gamma(idx) = 1; % Assign the randomly chosen states to 1
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

% Each lab recruits N_l patients, all with the same state gamma_l
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

rho_state = {}; % cell array where rho_state{i} is the kxk rho matrix for 
                % the ith state

for i=1:r
    rho_state{i} = randomCorrelationMatrix(sz);

%% Using simulated data (lab aggregates only), estimate correlation matrix