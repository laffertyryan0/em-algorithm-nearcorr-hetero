function out = vecL(mat)
    % Get the vectorized (flattened) part of the matrix strictly below
    % the diagonal. 
    out = mat(tril(true(size(mat)),-1));
end