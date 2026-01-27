function P = getMaskOrderingMatrix(mask_vector)
    [~,sortidx] = sort(mask_vector,'descend');
    numRow = length(sortidx);
    P = sparse(1:numRow,1:numRow,1);
    P = P(sortidx,:); %permutation corresponding to that sorting
end