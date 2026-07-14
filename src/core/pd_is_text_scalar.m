function tf = pd_is_text_scalar(value)
%ISTEXTSCALAR True for a character vector or scalar string when available.
% MATLAB R2016b compatible.

    tf = ischar(value);
    if tf
        return;
    end

    hasIsString = (exist('isstring', 'builtin') == 5) || ...
                  (exist('isstring', 'file') == 2);
    if hasIsString
        tf = isstring(value) && isscalar(value);
    end
end
