function config = pd_upgrade_app_config(config)
%PD_UPGRADE_APP_CONFIG Upgrade saved GUI settings to the current schema.

    if ~isstruct(config) || numel(config) ~= 1
        error('postdata:BadConfiguration', 'Configuration must be a scalar structure.');
    end
    if ~isfield(config, 'schemaVersion'), config.schemaVersion = 1; end
    if config.schemaVersion > 2
        error('postdata:NewerConfiguration', ...
            'Configuration schema %g is newer than this software.', config.schemaVersion);
    end
    config = addDefault(config, 'filePath', '');
    config = addDefault(config, 'task', 'chunk');
    config = addDefault(config, 'selectionMode', 'Index');
    config = addDefault(config, 'selectionValue', '1');
    config = addDefault(config, 'slurmPath', '');
    config = addDefault(config, 'progressMode', 'auto');
    config = addDefault(config, 'resultLevel', 'standard');

    computeDefault = pd_catalog_table_data(pd_option_catalog(config.task));
    plotDefault = pd_catalog_table_data(pd_plot_option_catalog());
    outputDefault = pd_catalog_table_data(pd_output_option_catalog());
    config.computeData = mergeData(computeDefault, getField(config, 'computeData', cell(0, 3)));
    config.plotData = mergeData(plotDefault, getField(config, 'plotData', cell(0, 3)));
    config.outputData = mergeData(outputDefault, getField(config, 'outputData', cell(0, 3)));
    config.schemaVersion = 2;
end

function output = mergeData(defaultData, savedData)
    output = defaultData;
    if ~iscell(savedData) || size(savedData, 2) < 2, return; end
    for i = 1:size(defaultData, 1)
        row = find(strcmp(defaultData{i, 1}, savedData(:, 1)), 1, 'first');
        if ~isempty(row), output{i, 2} = savedData{row, 2}; end
    end
end

function value = getField(input, name, defaultValue)
    if isfield(input, name), value = input.(name); else, value = defaultValue; end
end

function output = addDefault(output, name, value)
    if ~isfield(output, name), output.(name) = value; end
end
