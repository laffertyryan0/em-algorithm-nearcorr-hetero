dim = 5;
mat = patternedBlockCorrelation(dim, ...
                                .3, ...
                                20, ...
                                20,...
                                -.1);

% random walk but don't let condition number decrease
current_cond = 0;
while cond(mat) < 1200
    % randomly choose a entry to modify
    entryi = randi(dim);
    entryj= randi(dim);
    while entryi == entryj
        entryj= randi(dim);
    end
    delta = randn();
    mat_new = mat;
    mat_new(entryi,entryj) = mat_new(entryi,entryj) + delta;
    mat_new(entryj,entryi) = mat_new(entryi,entryj);
    new_cond = cond(mat_new);
    if new_cond>cond(mat) && min(eig(mat_new))>0 
        mat = mat_new;
    end
    disp(strcat("Condition Number: ",num2str(new_cond)));
end
disp(mat)