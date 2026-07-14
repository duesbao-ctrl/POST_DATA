function value = pd_parse_option_value(textValue, typeName)
%PD_PARSE_OPTION_VALUE Safely parse an editable option without EVAL/STR2NUM.

    if isnumeric(textValue) || islogical(textValue)
        raw = textValue;
    else
        raw = strtrim(pd_to_char(textValue));
    end

    switch lower(typeName)
        case 'numeric'
            if isnumeric(raw)
                if ~isreal(raw) || any(~isfinite(raw(:)))
                    error('postdata:BadNumericOption', ...
                        'Numeric options must contain only finite real values.');
                end
                value = raw;
                return;
            end
            if isempty(raw) || strcmp(raw, '[]')
                value = [];
                return;
            end
            cleaned = regexprep(raw, '[\[\],;]', ' ');
            tokens = regexp(strtrim(cleaned), '\s+', 'split');
            value = zeros(1, numel(tokens));
            for i = 1:numel(tokens)
                value(i) = str2double(tokens{i});
                if ~isfinite(value(i))
                    error('postdata:BadNumericOption', 'Invalid numeric value: %s', raw);
                end
            end

        case 'logical'
            if islogical(raw)
                if ~isscalar(raw)
                    error('postdata:BadLogicalOption', ...
                        'Logical options must be scalar.');
                end
                value = raw;
            elseif isnumeric(raw) && isscalar(raw) && isfinite(raw) && ...
                    any(raw == [0, 1])
                value = logical(raw);
            elseif any(strcmpi(raw, {'true','yes','on','1'}))
                value = true;
            elseif any(strcmpi(raw, {'false','no','off','0'}))
                value = false;
            else
                error('postdata:BadLogicalOption', 'Invalid logical value: %s', raw);
            end

        case 'celltext'
            if iscell(raw)
                value = raw;
                for i = 1:numel(value)
                    if ~pd_is_text_scalar(value{i})
                        error('postdata:BadCellTextOption', ...
                            'Cell-text options may contain text values only.');
                    end
                    value{i} = pd_to_char(value{i});
                end
            elseif isempty(raw) || strcmp(raw, '[]')
                value = {};
            else
                value = regexp(raw, '\s*,\s*', 'split');
                value = value(~cellfun(@isempty, value));
            end

        case 'text'
            if ~pd_is_text_scalar(raw)
                error('postdata:BadTextOption', 'Text option values must be text.');
            end
            raw = pd_to_char(raw);
            if ischar(raw) && strcmp(raw, '[]')
                value = '';
            else
                value = raw;
            end

        otherwise
            error('postdata:UnknownOptionType', 'Unknown option type: %s', typeName);
    end
end
