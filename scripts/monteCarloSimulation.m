%% Set path
addpath('../src')


%% Fix Parameters

num_metabolites = 50; %k
num_labs = 60; %L
average_fraction_missing_metabolites = 0.0;
num_mixture_components = 2; %r
mixing_probabilities = ones(1,num_mixture_components)/num_mixture_components;
num_subjects_per_lab = ones(num_labs,1)*100;  

rng(6); % Set seed for the whole simulation 
rho_state = {};
for j=1:num_mixture_components
    rho_state{j} = .2*ones(num_metabolites,num_metabolites) + ...
        .8*randomCorrelationMatrix(num_metabolites);
    if j==2
        rho_state{j} = eye(num_metabolites,num_metabolites);
    end
    if j==1
       rho_state{j} = patternedBlockCorrelation(num_metabolites, ...
                                            .7, ...
                                            20, ...
                                            20,...
                                            -.1);
    end
    assert(min(eig(rho_state{j}))>= 0, ...
        "Non-PSD matrix found for simulation rho_i.")
end
% Generate a random mask matrix for lab l
p = 1-average_fraction_missing_metabolites; %p = 0 means none are missing
for l=1:num_labs;
    subject_data_mask{l} = rand(1,num_metabolites) < p; % 1 if non-missing
end

%% Begin Monte Carlo Loop

final_rho_estimates = {};

MC_STEPS = 200;
for mc_step = 1:MC_STEPS
    
    [reported_spearman,...
              reported_spearman_mask,...
              n_samples,...
              r,...
              true_rho]...
              ...
                  =  simulateData(num_metabolites, ...
                                  num_labs, ...
                                  average_fraction_missing_metabolites, ...
                                  num_mixture_components, ...
                                  subject_data_mask,...
                                  rho_state,...
                                  mixing_probabilities,...
                                  num_subjects_per_lab,...
                                  []... % No seed for the single rep
                                  );
    
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
    
    % Calculate Pearson using Gaussian assumption
    reported_pearson = cellfun(@spearmanToPearson,reported_spearman, ...
                                            'UniformOutput',false);
    reported_pearson_vecL = cellfun(@vecL,reported_pearson,...
                                        'UniformOutput',false);
    
    
    % Calculate Fisher transformed correlation matrices
    reported_fisher = cellfun(@atanh,reported_pearson, ...
                                            'UniformOutput',false);
    reported_fisher_vecL = cellfun(@vecL,reported_fisher,...
                                        'UniformOutput',false);
    
    
    % For each l, find the permutation matrix P_l that puts missing entries 
    % below observed entries. Let [X Z]' = P_l * vecL(reported_pearson)
    P = {}; 
    X = {}; % Consists of all observed entries. 
            % Z here would be represented as all 0's due to the masking, but 
            % really Z means unobserved data. 
    
    for l=1:num_labs
        P{l} = getMaskOrderingMatrix(reported_spearman_mask_vecL{l});
        num_observed = sum(reported_spearman_mask_vecL{l});
        X_Z = P{l}*reported_fisher_vecL{l};
        X{l} = X_Z(1:num_observed);
    end
    
    
    MAX_EM_ITERATIONS = 30; % Outer loop
    MAX_GD_ITERATIONS = 1; % Inner PGD loop
    GD_TOLERANCE = 1;
    GD_LEARNING_RATE = 100*(.2/num_labs)/max(n_samples);
    %INIT_GDVARS_RANDLY = false;
    %NEARCORR_PROJ = false; % Do the correlation projection in the gd loop
    
    
    % Initialize EM parameters
    alpha_est = rand(1,r);
    alpha_est = alpha_est/sum(alpha_est); % Random initialization
    rho_est = cell(1,r);
    sigma_rho_est = cell(1,r);

    % Initialize EM parameters for NEARCORR_PROJ = false case
    alpha_est_no_proj = alpha_est;
    rho_est_no_proj = rho_est;
    sigma_rho_est_no_proj = sigma_rho_est;
    
    pearson_rho_est = cell(1,r); % Update this on every EM iteration
    
    for j = 1:r
        rho_est{j} = ...
            vecL(randomCorrelationMatrix(num_metabolites)); % Random initialization
        rho_est_no_proj{j} = ...
            rho_est{j}; % NCP = false
        sigma_rho_est{j} = ...
            speye(num_metabolites*(num_metabolites-1)/2);  
        sigma_rho_est_no_proj{j} = sigma_rho_est{j}; % NCP = false
    end
    
    w = 1./n_samples; % Lab-wise weighting factor for ...
                                     % variances (L vector)
    
    % Metrics to track for plotting. All should have prefix plotvar
    %plotvar_mse = {};  %rho mse
    %plotvar_bias = {}; %rho bias
    
    
    for em_iter=1:MAX_EM_ITERATIONS
        % disp("=========================================================")
        % fprintf("EM Iteration Number: %d\n",em_iter);
    
        % Log current estimate for mixing probabilities
        % fprintf("Current alpha estimate: ");
        % disp(alpha_est);
        % fprintf("\n");
    
        % Log intermediate spearman correlation matrix calculations
        % and compare with true pearson correlation
        % Log intermediate spearman correlation matrix calculations
        % and compare with true pearson correlation
    
            true_rho_fisher = cellfun(@atanh,true_rho, ... % the pearson
                                            'UniformOutput',false);
            rho_est_fisherinv = cellfun(@tanh,rho_est, ... % the pearson
                                            'UniformOutput',false);
            rho_est_fisherinv_no_proj = cellfun(@tanh,...
                                rho_est_no_proj, ... % the pearson
                                            'UniformOutput',false);

            estimated = {};
            estimated_no_proj = {};
            actual = {};    

            for j=1:r
                estimated{j} = vecLInverse(rho_est_fisherinv{1,j});
                estimated_no_proj{j} = vecLInverse(...
                                        rho_est_fisherinv_no_proj{1,j});
                actual{j} = true_rho{1,j};
            end
            
            order = inferComponentOrder(estimated,actual);
            order_no_proj = inferComponentOrder(estimated_no_proj,actual);
        
            % Append new values to plotvar metrics
            % for j=1:r
            %     estimated_j = rho_est_fisherinv{order(j)};
            %     actual_j = vecL(true_rho{j});
            %     plotvar_mse{j}(em_iter) = norm(estimated_j-actual_j,2);
            %     plotvar_bias{j}(em_iter) = mean(estimated_j-actual_j);
            % end
            % 
        % Show current alpha estimate
    
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
    
        [alpha_est,rho_est,sigma_rho_est] = argmaxQFisher(...
                                                    alpha_est,...
                                                    rho_est,...
                                                    sigma_rho_est, ...
                                                    X, ...
                                                    P, ...
                                                    w, ...
                                                    GD_LEARNING_RATE,...
                                                    MAX_GD_ITERATIONS, ...
                                                    GD_TOLERANCE, ...
                                                    false, ... %init rdly
                                                    true, ... %NEARCORR
                                                    em_iter);
        [alpha_est_no_proj,rho_est_no_proj,sigma_rho_est_no_proj] = ...
                                                    argmaxQFisher(...
                                                    alpha_est_no_proj,...
                                                    rho_est_no_proj,...
                                                 sigma_rho_est_no_proj, ...
                                                    X, ...
                                                    P, ...
                                                    w, ...
                                                    GD_LEARNING_RATE,...
                                                    MAX_GD_ITERATIONS, ...
                                                    GD_TOLERANCE, ...
                                                    false, ...%init rdly
                                                    false, ... %NEARCORR
                                                    em_iter);
    
    
       
    end
    final_rho_estimates{mc_step} = {};
    final_rho_estimates_no_proj{mc_step} = {};
    for j=1:r
        final_rho_estimates{mc_step}{j} = ...
            vecLInverse(rho_est_fisherinv{1,order(j)}); % last inferred order
        final_rho_estimates_no_proj{mc_step}{j} = ...
            vecLInverse(rho_est_fisherinv_no_proj...
                    {1,order_no_proj(j)}); % last inferred order
    end

    naive_rho = ...
        sum(cat(3,reported_pearson{:}),3)./...
        (sum(cat(3,reported_spearman_mask{:}),3)+eps); 
    while min(eig(naive_rho))<0
                naive_rho = naive_rho + .001*eye(size(naive_rho,1));
    end
            naive_rho = naive_rho / naive_rho(1,1);

    final_rho_estimates_naive{mc_step} = {};
    for j=1:r
        % Repeat this. It makes comparison easier
        final_rho_estimates_naive{mc_step}{j} = naive_rho;
    end
end

% Save run
save_filename = sprintf("../artifacts/mcruns/mcrun_%s",...
        string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
save(save_filename);


%% Load from saved workspace

% Set this to load a different workspace
LOAD_FILENAME = [];

% missing percent: .5
LOAD_FILENAME = "../artifacts/mcruns/mcrun_20260610_204433";

if ~isempty(LOAD_FILENAME)
    load(LOAD_FILENAME);
else
    load(save_filename);
end

%% Calculate metrics

% metrics:
% Bias
% MSE
% Median value of absolute difference (not outlier sensitive)
% Mean absolute deviation (from true value)
% Trimmed mean

% metric data structure
metrics = struct();
metrics.nearcorr = struct();
metrics.no_proj = struct();
metrics.naive = struct();

final_rho_estimates_all = {final_rho_estimates,...
                           final_rho_estimates_no_proj,...
                           final_rho_estimates_naive};

runtypes = fields(metrics);
for i=1:numel(runtypes)
    runtype = runtypes(i);
    runtype = runtype{:};
    metrics.(runtype) = struct();

    %%% TODO %%%
    metrics.(runtype).bias = {};
    metrics.(runtype).mse = {};
    metrics.(runtype).mae = {};
    metrics.(runtype).medae = {};
    metrics.(runtype).trim_mean = {};

    for j=1:r
        ests = [];
        for mc = 1:MC_STEPS
            ests = cat(3,ests,final_rho_estimates_all{i}{mc}{j});
        end
        metrics.(runtype).bias{j} = mean(ests-true_rho{1,j},3);
        metrics.(runtype).mse{j} = mean((ests-true_rho{1,j}).^2,3);
        metrics.(runtype).mae{j} = mean(abs(ests-true_rho{1,j}),3);
        metrics.(runtype).medae{j} = median(abs(ests-true_rho{1,j}),3);
        metrics.(runtype).trim_mean{j} = ...
                        trimmean(abs(ests-true_rho{1,j}),20,3);
    end
end


%% Create plots

num_rows = 8;
for j=1:r
    figure, 
    for m1=1:num_rows
       for m2=(m1+1):(num_rows+1)
            coef_ests = [];
            coef_ests_no_proj = [];
            coef_ests_naive = [];
            for step = 1:mc_step
                coef_act = true_rho{1,j}(m1,m2);
                coef_est = final_rho_estimates{step}{j};
                coef_est_no_proj = final_rho_estimates_no_proj{step}{j};
                coef_est_naive = final_rho_estimates_naive{step}{j};
                coef_ests = [coef_ests;...
                    coef_est(m1,m2)];
                coef_ests_no_proj = [coef_ests_no_proj;...
                    coef_est_no_proj(m1,m2)];
                coef_ests_naive = [coef_ests_naive;...
                    coef_est_naive(m1,m2)];
            end
            plotidx = (m1-1)*(num_rows+1) + m2;
            subplot(num_rows+1,num_rows+1,plotidx)
            hold on
           histogram(coef_ests-coef_act,'Normalization','probability', ...
               FaceColor='red');
            histogram(coef_ests_no_proj-coef_act,...
                'Normalization','probability',...
                FaceColor='blue');
            histogram(coef_ests_naive-coef_act,...
                'Normalization','probability',...
                FaceColor='green');
            hold off
            legend
        end           
    end
    sgtitle(strcat("Histogram of Correlation Estimates: ",...
                               " Component: ",...
                               string(j)))
end

for j=1:r
    figure,
    heatmap(metrics.nearcorr.mae{j}-metrics.no_proj.mae{j})
    title(strcat("Component ",string(j),...
        " -- MAE with proj - MAE without proj"));
    figure,
    heatmap(metrics.nearcorr.mse{j}-metrics.no_proj.mse{j})
    title(strcat("Component ",string(j),...
        " -- MSE with proj - MSE without proj"));
    figure,
    heatmap(metrics.nearcorr.trim_mean{j}-metrics.no_proj.trim_mean{j})
        title(strcat("Component ",string(j),...
            " -- 20% Trimmed MAE with proj - 20% Trimmed MAE without proj"));
end