function covR = estimateCovOfReportedFisher(reported_fisher, ...
                                                reported_spearman_mask, ...
                                                num_labs, ...
                                                n_samples)
    % Return the naive estimate for covariance of the correlations in 
    % vecL (flattened) form
    num_metabolites = length(reported_spearman_mask{1});
    num_pairs = num_metabolites*(num_metabolites-1)/2;
    rl_bar_weighted = sparse(zeros(num_pairs,1));
    N = sum(n_samples);
    denom = sparse(zeros(num_pairs,1));
    % Calculate total weighted mean
    for l=1:num_labs
        rl_bar_weighted = rl_bar_weighted + (n_samples(l))*...
                                        vecL(reported_fisher{l});
        denom = denom + ...
            n_samples(l)*vecL(reported_spearman_mask{l});
    end
    rl_bar_weighted = rl_bar_weighted ./ max(1,denom); % 0 if no data

    % Calculate weighted group level covariance (an unbiased estimate
    % of Sigma_rho = Sigma_rho{1} in one component case)
    covR = sparse(num_pairs,num_pairs);
    term_counts = sparse(zeros(num_pairs,1));
    for l=1:num_labs
        % rdiff = r_l - rbarbar 
        rdiff = vecL(reported_fisher{l}) - rl_bar_weighted;
        covR = covR + (n_samples(l)) * (rdiff*rdiff' ).*...
                        (vecL(reported_spearman_mask{l})*...
                        vecL(reported_spearman_mask{l})');
        % Only increase the count for an entry (i,j) if we have in lab L 
        % both i and j non-missing. that is if mask*mask' i,j = 1
        term_counts = term_counts + vecL(reported_spearman_mask{l})*...
            vecL(reported_spearman_mask{l})';
    end
    % We need to impute the variance if there's not enough data for any lab
    % 0 is ok for covariances, but for variances, if its 0, set it to 1
    %%% TODO %%%%
    covR(term_counts<=2 & logical(speye(num_pairs))) = 1;

    covR = covR ./ max(1,term_counts-1);

    % Ensure PSD

    while min(eig(covR))<0
        covR = speye(size(covR)).*covR + (1-speye(size(covR))).*covR * .5;
    end

    covR = sparse(covR);
