function options = pd_output_options_from_table(catalog, tableData)
%PD_OUTPUT_OPTIONS_FROM_TABLE Parse output settings from GUI table data.

    if size(tableData, 1) ~= numel(catalog)
        error('postdata:OutputOptionRowMismatch', 'Output table does not match its catalog.');
    end
    options = struct();
    for i = 1:numel(catalog)
        options.(catalog(i).Name) = pd_parse_option_value(tableData{i, 2}, catalog(i).Type);
    end
    if isempty(options.Directory)
        directoryRow = find(strcmp('Directory', {catalog.Name}), 1, 'first');
        options.Directory = pd_to_char(catalog(directoryRow).Default);
    end
    if ~(isscalar(options.DPI) && isfinite(options.DPI) && options.DPI > 0)
        error('postdata:BadDPI', 'DPI must be a positive scalar.');
    end
    dimensions = {'FigureWidthCm','FigureHeightCm'};
    for i = 1:numel(dimensions)
        value = options.(dimensions{i});
        if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
            error('postdata:BadFigureSize', ...
                '%s must be a positive finite scalar.', dimensions{i});
        end
    end
end
