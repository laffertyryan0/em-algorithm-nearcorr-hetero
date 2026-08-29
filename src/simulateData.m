function [reported_spearman,...
          reported_spearman_mask,...
          num_subjects_per_lab,...
          num_mixture_components,...
          true_rho]...
          ...
                      =  simulateData(num_metabolites, ...
                                      num_labs, ...
                                      average_fraction_missing_metabolites, ...
                                      num_mixture_components, ...
                                      subject_data_mask,...
                                      true_rho,...
                                      mixing_probabilities,...
                                      num_subjects_per_lab,...
                                      random_seed...
                                      )
    
    %set seed
    if ~isempty(random_seed)
        rng(random_seed);
    end
    
    tic
    
    % We have k metabolites
    k = num_metabolites;
    % We have L labs 
    L = num_labs;
    
    average_fraction_missing_metabolites = ...
          average_fraction_missing_metabolites;% approx what proportion of 
                                               % metabolites
                                               % are missing (1 = all missing
                                               % 0 = none missing)
                                               % sampled as iid bernoulli
    
    % Each lab has a state vector Gamma_l which is a one hot binary vector
    % Gamma_l,j = 1 iff lab is in state j. Let's say j = 1...r, maybe r=2 
    r = num_mixture_components;
    alpha = mixing_probabilities; % The alpha_j are the probabilities of a lab being 
                       % selected as state j
    
    
    Gamma = mnrnd(1,alpha,L); % multinomial dist with n=1 is categorical dist
                              % Gamma is Lxr with each row a One Hot Encoded
                              % lab state
    
    
    % Number of subjects in lab l is n_subjects(l)
    n_subjects = num_subjects_per_lab;
    
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
        rho_state{j} = true_rho{j};
    end
    clear j
    
    % Generate subject level data for each lab
    subject_data = {}; % subject_data{i} is design matrix for 
                             % ith lab
    % Define a mask for each lab: 1 for non-missing, 0 for missing
    % subject_data_mask{l} is a binary vector of length n_subjects(l)
    % Lab l only "really" collects data for unmasked (mask=1) metabolites
    % 1 = non-missing, 0 = missing
    
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
    
    
        covariate_mask = subject_data_mask{l};
        
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

    reported_spearman = reported_spearman;
    reported_spearman_mask = reported_spearman_mask;
    num_subjects_per_lab = num_subjects_per_lab;
    num_mixture_components = num_mixture_components; %r
    true_rho = rho_state;