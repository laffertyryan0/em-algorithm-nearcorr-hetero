%% Set path
addpath('../src')

%% Loop over missingness percentages

metrics_by_missingness = {};

missing_iters = 30; %.0 to .9 = 1 to 10
max_mi = 1;
missingnesses = ((1:missing_iters)-1)*(max_mi/missing_iters);
for missing_iter = 1:missing_iters
    
    %% Fix Parameters
    
    num_metabolites = 50; %k
    num_labs = 60; %L
    average_fraction_missing_metabolites = missingnesses(missing_iter);
    num_mixture_components = 1; %r
    mixing_probabilities = ones(1,num_mixture_components)/num_mixture_components;
    num_subjects_per_lab = ones(num_labs,1)*1000;
    
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

%%%DEBUG%%%%
    % % Ensure at least one lab per pair
    % for i1 = 1:num_metabolites
    %     for i2 = 1:num_metabolites
    %         in_at_least_one = false;
    %         for l=1:num_labs
    %             if subject_data_mask{l}(i1)==1 && ...
    %                     subject_data_mask{l}(i2)==1
    %                 in_at_least_one = true;
    %             end
    %         end
    %         if ~in_at_least_one
    %             % Randomly choose a lab
    %             forceL = randi(num_labs);
    %             % Force forceL to have both
    %             subject_data_mask{forceL}(i1) = 1;
    %             subject_data_mask{forceL}(i2) = 1;
    %         end
    %     end
    % end
%%%%%%%%%%%

    % Ensure at least two non-missing in each lab
    for l=1:num_labs
        if sum(subject_data_mask{l}) < 2
            for q=randsample(num_metabolites,2)
                subject_data_mask{l}(q) = 1;
            end
        end
    end
    
    %% Begin Monte Carlo Loop
    
    final_rho_estimates = {};

    %%%%%%%%%%% TEMPORARY -- FIX THIS %%%%%%%%%%
    CIs = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MC_STEPS = 20;
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
        
        
        MAX_EM_ITERATIONS = 500; % Outer loop
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

        % 
        % %%%%%%%%%%%%%%%%%%%%%%% TEMPORARY -- FIX THIS %%%%%%%%%%%%%%%%%%%
        % 
        %     FI_est = louisFisherEst(...
        %                                         alpha_est,...
        %                                         rho_est,...
        %                                         sigma_rho_est, ...
        %                                         X, ...
        %                                         P, ...
        %                                         w, ...
        %                                         GD_LEARNING_RATE,...
        %                                         MAX_GD_ITERATIONS, ...
        %                                         GD_TOLERANCE, ...
        %                                         false, ...
        %                                         0, ...
        %                                         em_iter);
        %     rho_FI_est = FI_est(max(1,r-1):end,max(1,r-1):end); % only rho part
        %     a = .05; 
        %     std_err = vecLInverse(sqrt(abs(diag(inv(rho_FI_est)))));
        %     %disp(std_err');
        %     j=1;
        %     CI_upper = final_rho_estimates{mc_step}{j} + norminv(1-a/2)*std_err;
        %     CI_lower = final_rho_estimates{mc_step}{j} - norminv(1-a/2)*std_err;
        %     CI = {CI_lower CI_upper};
        %     %disp(CI)
        %     CIs{missing_iter}{mc_step} = CI;
        % 
        % 
        % %%%%%%%%%%%%%%%%%%%%%%% TEMPORARY -- FIX THIS %%%%%%%%%%%%%%%%%%%
        % 
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
    metrics_by_missingness{missing_iter} = metrics;
end


%% Save run

save_filename = sprintf("../artifacts/mc_ml_runs/mc_ml_run_%s",...
        string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
save(save_filename);


%% Load from saved workspace

% Set this to load a different workspace
LOAD_FILENAME = [];

%LOAD_FILENAME = "../artifacts/mc_ml_runs/mc_ml_run_20260627_160242.mat";
%LOAD_FILENAME = "../artifacts/mc_ml_runs/mc_ml_run_20260812_002107.mat" %30
%missingnesses (one component)
%LOAD_FILENAME = "../artifacts/mc_ml_runs/mc_ml_run_20260812_123859.mat" %30
%missingnesses (two component)
%LOAD_FILENAME = "../artifacts/mc_ml_runs/mc_ml_run_20260818_232835.mat"
% ^ test run with CIs (not fully complete code)

LOAD_FILENAME = "../artifacts/mc_ml_runs/mc_ml_run_20260821_111741.mat"

if ~isempty(LOAD_FILENAME)
    disp("Using saved file.")
    load(LOAD_FILENAME);
else
    load(save_filename);
end

%% Create plots

% First a plot of the (a,b)th component MSE over probability

a = 3;
b = 2;

mse_by_mprob_nearcorr = [];
mse_by_mprob_no_proj = [];
mse_by_mprob_naive = [];

for j=1:r
    figure,
    
    for i=1:missing_iters
        mse_by_mprob_nearcorr(i) = ...
            mean(metrics_by_missingness{i}.nearcorr.mae{j},'all');
        mse_by_mprob_no_proj(i) = ...
            mean(metrics_by_missingness{i}.no_proj.mae{j},'all');
        mse_by_mprob_naive(i) = ...
            mean(metrics_by_missingness{i}.naive.mae{j},'all');

        % 
        % mse_by_mprob_nearcorr(i) = ...
        %     metrics_by_missingness{i}.nearcorr.medae{j}(a,b);
        % mse_by_mprob_no_proj(i) = ...
        %     metrics_by_missingness{i}.no_proj.medae{j}(a,b);
        % mse_by_mprob_naive(i) = ...
        %     metrics_by_missingness{i}.naive.medae{j}(a,b);
    end
    
    hold on;
    plot(missingnesses,mse_by_mprob_nearcorr,'Color','blue','DisplayName','EM with projection');
    plot(missingnesses,mse_by_mprob_no_proj,'Color','red','DisplayName','EM w/o projection');
    % Only plot the naive if there are more than one component
    % Otherwise its too confusing
    if r<=1
        plot(missingnesses,mse_by_mprob_naive,'Color','green','DisplayName','Naive');
    end
    title(strcat("Component: ",string(j)))
    legend('Location','best')
    hold off;
end

% Then create a plot where all the components are on the same plot

figure,
hold on
grid on
%title("ALL on same plot")

for j=1:r
    figure,
    for a=1:(num_metabolites-1)
        for b=(a+1):num_metabolites
            
            mse_by_mprob_nearcorr = [];
            mse_by_mprob_no_proj = [];
            mse_by_mprob_naive = [];
            
                
                for i=1:missing_iters
                    % mse_by_mprob_nearcorr(i) = ...
                    %     mean(metrics_by_missingness{i}.nearcorr.medae{j},'all');
                    % mse_by_mprob_no_proj(i) = ...
                    %     mean(metrics_by_missingness{i}.no_proj.medae{j},'all');
                    % mse_by_mprob_naive(i) = ...
                    %     mean(metrics_by_missingness{i}.naive.medae{j},'all');
            
            
                    mse_by_mprob_nearcorr(i) = ...
                        metrics_by_missingness{i}.nearcorr.bias{j}(a,b);
                    mse_by_mprob_no_proj(i) = ...
                        metrics_by_missingness{i}.no_proj.bias{j}(a,b);
                    mse_by_mprob_naive(i) = ...
                        metrics_by_missingness{i}.naive.bias{j}(a,b);
                end
                hold on
                p=plot(missingnesses*100,mse_by_mprob_nearcorr,'Color','blue',...
                    'LineWidth',.02);
                p.Color(4) = .05;
                p2=plot(missingnesses*100,mse_by_mprob_no_proj,'Color','red',...
                    'LineWidth',.02);
                p2.Color(4) = .05;           
                if r<=1
                    p3=plot(missingnesses*100,mse_by_mprob_naive,'Color','green',...
                       'LineWidth',.02);
                    p3.Color(4) = .05;      
                end
               % p.Color(4) = .3;

                hold off
                % p.Color(4) = .5;
               % plot(missingnesses,mse_by_mprob_naive,'Color','green',...
               %     'LineWidth',.2);
               % p.Color(4) = .3;
                % xlim([0,.6]);
                % ylim([-.05,.05]);
    

                xlim([0,100]);
                ylim([-.5,.5]);
                ylabel("Estimate - True")
                xlabel("Average Percent Variables Missing ")
               % title("Monte-Carlo Bias Estimates for all Pairwise Correlations")
                title(strcat("Component: ",string(j)))
            
        end
    end
hold off
set(gcf,'Color','white')
set(gcf, 'Position', [100, 100, 800, 400])
set(gcf, 'InvertHardcopy', 'off'); 
set(gcf, 'PaperPositionMode', 'auto');
end

% mc_ml_run_20260812_002107.mat -- for 1 component with 30 missingnesses,
% 20 mc runs per iteration, 500 EM steps

%% TEMPORARY -- remove this -- this is debug code for CIs
% coverage prob = percent of cis that contain the true rho
     %   i = 48; j = 4; % look at (i,j) correlation component
% for missing_iter = 3:missing_iters
%     for i=1:num_metabolites
%         for j=1:num_metabolites
%             total = 0;
%             for mc_step = 1:MC_STEPS
%                 if CIs{missing_iter}{mc_step}{1}(i,j) <= true_rho{1}(i,j) && ...
%                         true_rho{1}(i,j) <= CIs{missing_iter}{mc_step}{2}(i,j)
%                     total = total + 1;
%                 end
%             end
%             cov_prob(i,j) = total/length(CIs{missing_iter});
%         end
%     end
%     figure,
%     title(strcat("Missingness: ",string(missingnesses(missing_iter))))
%     heatmap(cov_prob)
% end