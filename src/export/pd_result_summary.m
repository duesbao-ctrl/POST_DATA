function data = pd_result_summary(result)
%PD_RESULT_SUMMARY Flatten scalar result metadata for tables and CSV output.

    rows = cell(0, 2);
    rows = appendStruct(rows, result, '', 0);
    data = rows;
end

function rows = appendStruct(rows, value, prefix, depth)
    if depth > 6 || ~isstruct(value) || numel(value) ~= 1
        return;
    end
    fields = fieldnames(value);
    for i = 1:numel(fields)
        name = fields{i};
        itemValue = value.(name);
        if isempty(prefix)
            path = name;
        else
            path = [prefix, '.', name];
        end
        if isstruct(itemValue)
            rows = appendStruct(rows, itemValue, path, depth + 1);
        elseif islogical(itemValue) && isscalar(itemValue)
            rows(end + 1, :) = {path, pd_format_option_value(itemValue)}; %#ok<AGROW>
        elseif isnumeric(itemValue) && isscalar(itemValue) && ~ishghandle(itemValue)
            rows(end + 1, :) = {path, pd_format_option_value(itemValue)}; %#ok<AGROW>
        elseif ischar(itemValue) && size(itemValue, 1) == 1 && numel(itemValue) <= 240
            rows(end + 1, :) = {path, itemValue}; %#ok<AGROW>
        end
    end
end
