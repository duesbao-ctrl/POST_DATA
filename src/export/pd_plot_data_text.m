function text = pd_plot_data_text(plotData, delimiter, rowIndices, includeHeader)
%PD_PLOT_DATA_TEXT Serialize plot-view numeric data as delimited text.

    if nargin < 2 || isempty(delimiter)
        delimiter = ',';
    end
    if nargin < 3 || isempty(rowIndices)
        rowIndices = 1:size(plotData.Values, 1);
    end
    if nargin < 4 || isempty(includeHeader)
        includeHeader = true;
    end
    if ~(islogical(includeHeader) && isscalar(includeHeader)) && ...
            ~(isnumeric(includeHeader) && isscalar(includeHeader) && ...
            isfinite(includeHeader) && any(includeHeader == [0 1]))
        error('postdata:plotDataText:BadHeaderFlag', ...
            'includeHeader must be a logical scalar.');
    end
    includeHeader = logical(includeHeader);
    validateContract(plotData);
    rowIndices = rowIndices(:).';
    if any(rowIndices < 1) || any(rowIndices > size(plotData.Values, 1)) || ...
            any(rowIndices ~= round(rowIndices))
        error('postdata:plotDataText:BadRows', ...
            'Row indices must select existing integer rows.');
    end

    lineCount = numel(rowIndices) + double(includeHeader);
    lines = cell(lineCount, 1);
    offset = 0;
    if includeHeader
        header = cell(size(plotData.ColumnNames));
        for i = 1:numel(header)
            header{i} = escapeText(plotData.ColumnNames{i}, delimiter);
        end
        lines{1} = strjoin(header, delimiter);
        offset = 1;
    end
    for i = 1:numel(rowIndices)
        values = plotData.Values(rowIndices(i), :);
        cells = cell(1, numel(values));
        for j = 1:numel(values)
            cells{j} = sprintf('%.16g', values(j));
        end
        lines{i + offset} = strjoin(cells, delimiter);
    end
    text = strjoin(lines, sprintf('\n'));
    text = [text, sprintf('\n')];
end

function value = escapeText(value, delimiter)
    value = pd_to_char(value);
    value = strrep(value, sprintf('\r'), ' ');
    value = strrep(value, sprintf('\n'), ' ');
    quote = ~isempty(strfind(value, delimiter)) || ... %#ok<STREMP>
        ~isempty(strfind(value, '"')); %#ok<STREMP>
    if quote
        value = strrep(value, '"', '""');
        value = ['"', value, '"'];
    end
end

function validateContract(plotData)
    required = {'ColumnNames','Values'};
    if ~isstruct(plotData) || ~all(isfield(plotData, required))
        error('postdata:plotDataText:BadContract', ...
            'Plot data must contain ColumnNames and Values.');
    end
    if ~isnumeric(plotData.Values) || ndims(plotData.Values) ~= 2
        error('postdata:plotDataText:BadValues', ...
            'Plot data Values must be a numeric matrix.');
    end
    if numel(plotData.ColumnNames) ~= size(plotData.Values, 2)
        error('postdata:plotDataText:BadColumns', ...
            'ColumnNames must match the number of data columns.');
    end
end
