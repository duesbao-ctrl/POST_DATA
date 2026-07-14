function mode = pd_distribution_normalize_empty_bin_mode(value)
%PD_DISTRIBUTION_NORMALIZE_EMPTY_BIN_MODE Normalize empty-bin policy.

    mode = lower(strtrim(toChar(value)));
    switch mode
        case {'zero', '0'}
            mode = 'zero';
        case {'nan', 'na'}
            mode = 'nan';
        case {'remove', 'delete', 'drop'}
            mode = 'remove';
    end
end

function out = toChar(value)
    if ischar(value)
        out = value;
    else
        out = char(value);
    end
end
