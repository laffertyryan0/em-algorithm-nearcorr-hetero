function out = randomCorrelationMatrix(sz, ...
                                       off_diag_const,...
                                       corr_density_param,... % in [0,1]
                                       random_seed)

    if ~isempty(random_seed)
        rng(random_seed)
    end

    vec = rand(sz,1)<corr_density_param;
    mat = off_diag_const * vec * vec';
    mat(logical(eye(length(mat)))) = 1;
    
    out = mat
end
    