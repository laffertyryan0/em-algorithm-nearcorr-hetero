function out = vecL(mat)
    out = mat(tril(true(size(mat)),-1));
end