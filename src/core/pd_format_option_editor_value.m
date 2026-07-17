function textValue = pd_format_option_editor_value(value, typeName)
%PD_FORMAT_OPTION_EDITOR_VALUE Format a value for an editable option cell.
% Empty values are blank, numeric scalars have no brackets, and numeric
% vectors retain brackets so their shape remains visually explicit.

    if nargin < 2
        typeName = '';
    end
    typeName = lower(strtrim(pd_to_char(typeName)));

    if (ischar(value) || pd_is_text_scalar(value)) && ...
            any(strcmp(typeName, {'numeric','logical','celltext'}))
        try
            value = pd_parse_option_value(value, typeName);
        catch
            textValue = pd_to_char(value);
            return;
        end
    end

    if isempty(value)
        textValue = '';
    else
        textValue = pd_format_option_value(value);
    end
end
