function out = vecLInverse(v)
    % Get the symmetric matrix formed by placing v in the lower triangular
    % part and filling the diagonal with ones.
    d = floor((1+sqrt(1+8*length(v)))/2); % Pos. solution of x = d(d-1)/2
    out = zeros(d,d);
    out(tril(true(d,d),-1))=v;
    out = out + out' + eye(d,d);
end