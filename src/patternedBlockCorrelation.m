function out = patternedBlockCorrelation(sz, ...,
                                         off_diag_const,...
                                         pattern_block_size, ...
                                         pattern_block_overlap,...
                                         outside_pattern)
    out = zeros(sz)+outside_pattern;
    for i=1:pattern_block_overlap:sz
        if (i+pattern_block_size)<=sz
            out(i:(i+pattern_block_size-1),i:(i+pattern_block_size-1)) = ...
                off_diag_const;
        else
            out(i:end,i:end) = off_diag_const;
        end

    out(logical(eye(length(out)))) = 1;
end
    