function value = pd_to_char(value)
%TOCHAR Convert a scalar string to char without requiring modern syntax.
% MATLAB R2016b compatible.

    hasIsString = (exist('isstring', 'builtin') == 5) || ...
                  (exist('isstring', 'file') == 2);
    if hasIsString && isstring(value)
        value = char(value);
    end
end
