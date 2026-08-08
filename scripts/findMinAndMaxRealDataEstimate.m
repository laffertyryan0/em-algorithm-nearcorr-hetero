addpath("../src")

est = load('../artifacts/rho_est_fisherinv_nearcorrtrue.mat');
mat = vecLInverse(rho_est_fisherinv{1});
mat(logical(eye(size(mat,1)))) = 0;
[~,idx] = max(mat,[],'all','linear');
[r,c] = ind2sub(size(mat),idx);
disp(strcat("max Row: ",string(r)))
disp(strcat("max Col: ",string(c)))

est = load('../artifacts/rho_est_fisherinv_nearcorrtrue.mat');
mat = vecLInverse(rho_est_fisherinv{1});
mat(logical(eye(size(mat,1)))) = 0;
[~,idx] = min(mat,[],'all','linear');
[r,c] = ind2sub(size(mat),idx);
disp(strcat("min Row: ",string(r)))
disp(strcat("min Col: ",string(c)))

est = load('../artifacts/rho_est_fisherinv_nearcorrtrue.mat');
mat = vecLInverse(rho_est_fisherinv{1});
mat(logical(eye(size(mat,1)))) = Inf;
[~,idx] = min(abs(mat),[],'all','linear');
[r,c] = ind2sub(size(mat),idx);
disp(strcat("minabs Row: ",string(r)))
disp(strcat("minabs Col: ",string(c)))

