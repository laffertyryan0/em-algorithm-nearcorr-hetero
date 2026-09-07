min_eigs = [];
norm_delta = []; % improvement in norm when you do nearcorr
norm_delta_nomissingest = [];
mat_eig_viol = [];
for trial = 1:1000
    mat = rand(5,5).*(1-eye(5,5)) + eye(5,5);
    mat = nearcorr((mat + mat')/2);

    if min(eig(mat))<0
        mat_eig_viol(trial) = 1;
    else
        mat_eig_viol(trial) = 0;
    end
    
    num_samples = 500;
    
    data = zeros(num_samples,5,2);
    for l=1:3
        data(:,:,l) = mvnrnd(zeros(5,1),mat,num_samples);
    end
        
        corrnomask = {};
        for i=1:3
            corrnomask{i} = corr(data(:,:,i));
        end
        
        nomissingest = mean(cat(3,corrnomask{:}),3);
    
    
        data(:,1,1) = nan;
    
        data(:,5,2) = nan;
        data(:,3,3) = nan;
    
    % How can we make up an estimate by using missing entries in such a way
    % that the naive correlation estimate is not pd
    
    corrs = {};
    
    corrs{1} = corr(data(:,:,1));
    corrs{2} = corr(data(:,:,2));
    corrs{3} = corr(data(:,:,3));
    
    corrsarr = cat(3,corrs{:});
    
    corr_est = zeros(5,5);
    counts = zeros(5,5);
    for i=1:5
        for j=1:5
            counts(i,j) = 0;
            for l=1:3
                delta = corrsarr(i,j,l);
                if ~isnan(delta)
                    corr_est(i,j) = corr_est(i,j) + delta;
                    counts(i,j) = counts(i,j) + 1;
                end
            end
            corr_est(i,j) = corr_est(i,j)/counts(i,j);
        end
    end
    
    % disp(min(eig(corr_est)))
    % disp("norm(mat-corr_est)")
    % disp(norm(mat-corr_est))
    % disp("norm(mat-nearcorr(corr_est))")
    % disp(norm(mat-nearcorr(corr_est)))
    min_eigs(trial) = min(eig(corr_est));
    norm_delta(trial) = -norm(mat-corr_est,'fro') + norm(mat-nearcorr(corr_est),'fro');
    norm_delta_nomissingest(trial) = -norm(nomissingest-corr_est,'fro') + ...
                    norm(nomissingest-nearcorr(corr_est),'fro');
    


end

figure,
scatter(min_eigs,norm_delta,1);
xlabel("Minimum Eigenvalue of Unprojected Estimate")
ylabel("Change in Frobenius distance to true correlation")

figure,
hist(min_eigs)

% figure,
% scatter(min_eigs,norm_delta_nomissingest,1);
% xlabel("Minimum Eigenvalue of Unprojected Estimate")
% ylabel("Change in Frobenius distance to unmasked estimate")