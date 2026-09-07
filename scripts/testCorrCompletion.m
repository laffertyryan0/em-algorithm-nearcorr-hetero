completions = []
min_eigs = []
log_dets = []
mat = rand(5,5).*(1-eye(5,5)) + eye(5,5);
mat = nearcorr((mat + mat')/2);
for i=1:100
    a  = mat;%final_rho_estimates_no_proj{1, 1}{1, 1}; % modified
    b = a; % original
    a(1,2) = rand(1)*2-1;
    a(2,1) = a(1,2);
    % disp("EIG (A)")
    % disp(eig(a))
    w = ones(5,5);
    w(1,2) = 0;
    w(2,1) = 0;
    c = nearcorr(a,"Weights",w); % re estimated
    % disp("||B-C||")
    % disp(norm(b-c));
    % disp("||B-A||")
    % disp(norm(b-a));
    % disp("C-B")
    % disp(c-b);
    completions(i) = c(1,2);
    min_eigs(i) = min(eig(c));
    log_dets(i) = log(det(c));
 end

 figure,
 histogram(completions,30)
 figure,
 scatter(min_eigs,completions)
 xlabel("Min eigs")
figure,
 scatter(log_dets,completions)
 xlabel("Log dets")