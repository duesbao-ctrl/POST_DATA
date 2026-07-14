function options = pd_plot_options_from_table(catalog, tableData)
%PD_PLOT_OPTIONS_FROM_TABLE Parse plot settings from GUI table data.

    if size(tableData, 1) ~= numel(catalog)
        error('postdata:PlotOptionRowMismatch', 'Plot table does not match its catalog.');
    end
    options = struct();
    for i = 1:numel(catalog)
        value = pd_parse_option_value(tableData{i, 2}, catalog(i).Type);
        if strcmp(catalog(i).Name, 'UpdateMode')
            modeAlias = lower(strtrim(pd_to_char(value)));
            if any(strcmp(modeAlias, {'overwrite'})), value = 'replace'; end
            if any(strcmp(modeAlias, {'add','append'})), value = 'overlay'; end
            if ~any(strcmp(value, {'replace','overlay'}))
                error('postdata:BadPlotUpdateMode', ...
                    'UpdateMode must be replace or overlay.');
            end
        end
        if ~isempty(catalog(i).AllowedValues)
            textValue = pd_to_char(value);
            index = find(strcmpi(textValue, catalog(i).AllowedValues), 1, 'first');
            if isempty(index)
                error('postdata:BadPlotChoice', ...
                    '%s must be one of: %s.', catalog(i).Name, ...
                    joinChoices(catalog(i).AllowedValues));
            end
            value = catalog(i).AllowedValues{index};
        end
        options.(catalog(i).Name) = value;
    end
    if ~isfield(options, 'UpdateMode')
        options.UpdateMode = 'replace';
    end
    mode = lower(strtrim(pd_to_char(options.UpdateMode)));
    switch mode
        case {'replace', 'overwrite'}
            options.UpdateMode = 'replace';
        case {'overlay', 'add', 'append'}
            options.UpdateMode = 'overlay';
        otherwise
            error('postdata:BadPlotUpdateMode', ...
                'UpdateMode must be replace or overlay.');
    end
    positive = {'FontSize','TitleFontSize','AxisLineWidth','LineWidth', ...
        'MarkerSize','ContourLevels'};
    for i = 1:numel(positive)
        value = options.(positive{i});
        if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
            error('postdata:BadPlotNumericOption', ...
                '%s must be a positive finite scalar.', positive{i});
        end
    end
    if abs(options.ContourLevels - round(options.ContourLevels)) > eps || ...
            options.ContourLevels < 2
        error('postdata:BadContourLevels', ...
            'ContourLevels must be an integer greater than or equal to 2.');
    end
    if ~(isscalar(options.GridAlpha) && isfinite(options.GridAlpha) && ...
            options.GridAlpha >= 0 && options.GridAlpha <= 1)
        error('postdata:BadGridAlpha', 'GridAlpha must be between 0 and 1.');
    end
    ranges = {'ColorLimits','XLim','YLim'};
    for i = 1:numel(ranges)
        value = options.(ranges{i});
        if ~isempty(value) && ~(isnumeric(value) && numel(value) == 2 && ...
                all(isfinite(value)) && value(2) > value(1))
            error('postdata:BadPlotRange', ...
                '%s must be empty or [min max].', ranges{i});
        end
    end
    colors = options.CustomColorOrder;
    if ~isempty(colors) && ~(isnumeric(colors) && ismatrix(colors) && ...
            size(colors, 2) == 3 && all(isfinite(colors(:))) && ...
            all(colors(:) >= 0) && all(colors(:) <= 1))
        error('postdata:BadCustomColorOrder', ...
            'CustomColorOrder must be empty or an N-by-3 RGB matrix from 0 to 1.');
    end
end

function text = joinChoices(values)
    text = values{1};
    for i = 2:numel(values), text = [text, ', ', values{i}]; end %#ok<AGROW>
end
