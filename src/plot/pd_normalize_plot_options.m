function options = pd_normalize_plot_options(input)
%PD_NORMALIZE_PLOT_OPTIONS Validate partial plot settings and fill defaults.

    catalog = pd_plot_option_catalog();
    data = pd_catalog_table_data(catalog);
    if nargin < 1 || isempty(input)
        options = pd_plot_options_from_table(catalog, data);
        return;
    end
    if isstruct(input) && isscalar(input)
        names = fieldnames(input);
        for i = 1:numel(names)
            row = find(strcmpi(names{i}, {catalog.Name}), 1, 'first');
            if isempty(row)
                error('postdata:UnknownPlotOption', ...
                    'Unknown plot option: %s.', names{i});
            end
            data{row, 2} = input.(names{i});
        end
    elseif iscell(input) && mod(numel(input), 2) == 0
        for i = 1:2:numel(input)
            if ~pd_is_text_scalar(input{i})
                error('postdata:BadPlotOptionName', ...
                    'Plot option names must be text.');
            end
            name = pd_to_char(input{i});
            row = find(strcmpi(name, {catalog.Name}), 1, 'first');
            if isempty(row)
                error('postdata:UnknownPlotOption', ...
                    'Unknown plot option: %s.', name);
            end
            data{row, 2} = input{i + 1};
        end
    else
        error('postdata:BadPlotStyleOptions', ...
            'Plot style options must be a scalar struct or name-value cell array.');
    end
    options = pd_plot_options_from_table(catalog, data);
end
