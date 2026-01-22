function out = randomCorrelationMatrix(sz)
    % size sz
    % out = a random correlation matrix of size sz x sz
    A = randn(sz, sz);
    S = A' * A;
    d = sqrt(diag(S)); 
    out = S ./ (d * d'); 
    out(1:sz+1:end) = 1;