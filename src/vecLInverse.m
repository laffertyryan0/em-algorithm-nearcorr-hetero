function out = vecLInverse(v)
    d = floor((1+sqrt(1+8*length(v)))/2); % Pos. solution of x = d(d-1)/2
    out = zeros(d,d);
    out(tril(true(d,d),-1))=v;
    out = out + out' + eye(d,d);
end