function textValue = pd_format_option_value(value)
%PD_FORMAT_OPTION_VALUE Convert an option value to editable table text.

    if isempty(value)
        textValue = '[]';
    elseif islogical(value)
        if value
            textValue = 'true';
        else
            textValue = 'false';
        end
    elseif isnumeric(value)
        if isscalar(value)
            textValue = num2str(value, 16);
        else
            parts = arrayfun(@(x) num2str(x, 16), value(:).', 'UniformOutput', false);
            textValue = ['[', strjoin(parts, ' '), ']'];
        end
    elseif iscell(value)
        textValue = strjoin(value, ',');
        if isempty(textValue)
            textValue = '[]';
        end
    else
        textValue = pd_to_char(value);
    end
end
