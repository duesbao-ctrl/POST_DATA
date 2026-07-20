function [values, enabled] = pd_validate_catalog_values(catalog, tableData)
%PD_VALIDATE_CATALOG_VALUES Parse and semantically validate active options.

    if size(tableData, 1) ~= numel(catalog)
        error('postdata:OptionRowMismatch', 'Option table does not match the active catalog.');
    end
    enabled = pd_catalog_enabled_mask(catalog, tableData);
    values = cell(numel(catalog), 1);
    for i = 1:numel(catalog)
        if ~enabled(i), continue; end
        try
            value = pd_parse_option_value(tableData{i, 2}, catalog(i).Type);
            validateOne(catalog(i), value);
            values{i} = value;
        catch err
            error('postdata:InvalidOption', '%s: %s', catalog(i).Name, err.message);
        end
    end
end

function validateOne(meta, value)
    if meta.Required && isempty(value)
        error('A value is required.');
    end
    if ~isempty(meta.AllowedValues) && ~isempty(value)
        if isnumeric(value)
            allowed = cellfun(@str2double, meta.AllowedValues);
            valid = isscalar(value) && any(value == allowed);
        else
            valid = ischar(value) && any(strcmpi(value, meta.AllowedValues));
        end
        if ~valid
            error('Allowed values: %s.', strjoin(meta.AllowedValues, ', '));
        end
    end
    if isnumeric(value) && ~isempty(value)
        if ~isempty(meta.Minimum) && meta.ExclusiveMinimum && any(value <= meta.Minimum)
            error('Value must be greater than %g.', meta.Minimum);
        elseif ~isempty(meta.Minimum) && any(value < meta.Minimum)
            error('Value must be at least %g.', meta.Minimum);
        end
        if ~isempty(meta.Maximum) && any(value > meta.Maximum)
            error('Value must not exceed %g.', meta.Maximum);
        end
        if meta.IntegerOnly && any(abs(value - round(value)) > 1e-12)
            error('Value must be an integer.');
        end
        if ~isempty(meta.VectorLength) && numel(value) ~= meta.VectorLength
            error('Value must contain %d numbers.', meta.VectorLength);
        end
        if isequal(meta.VectorLength, 2) && value(2) < value(1)
            error('Range maximum must be greater than or equal to minimum.');
        end
    end
end
