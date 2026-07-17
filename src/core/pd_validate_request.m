function request = pd_validate_request(request)
%VALIDATEREQUEST Validate and normalize a versioned POST_DATA request.
% MATLAB R2016b compatible.

    if ~isstruct(request) || ~isscalar(request)
        error('postdata:validateRequest:BadRequest', 'Request must be a scalar structure.');
    end
    required = {'schemaVersion', 'analysisType', 'baseDir', 'filePath', ...
        'selection', 'progressMode', 'makePlots', 'chunk', 'analysisOptions'};
    for i = 1:numel(required)
        if ~isfield(request, required{i})
            error('postdata:validateRequest:MissingField', ...
                'Request is missing field "%s".', required{i});
        end
    end
    if ~(isnumeric(request.schemaVersion) && isscalar(request.schemaVersion) && ...
            isfinite(request.schemaVersion) && request.schemaVersion == 1)
        error('postdata:validateRequest:SchemaVersion', ...
            'Unsupported request schema version. Expected numeric value 1.');
    end

    if ~pd_is_text_scalar(request.analysisType)
        error('postdata:validateRequest:BadType', ...
            'analysisType must be a text scalar.');
    end
    request.analysisType = pd_normalize_analysis_type(request.analysisType);
    validTypes = {'chunk', 'cluster', 'vx', 'massx', 'network2d'};
    if ~any(strcmp(request.analysisType, validTypes))
        error('postdata:validateRequest:BadType', 'Unsupported analysis type.');
    end
    request.baseDir = pd_to_char(request.baseDir);
    request.filePath = pd_to_char(request.filePath);
    request.progressMode = pd_to_char(request.progressMode);
    if ~pd_is_text_scalar(request.baseDir) || ~pd_is_text_scalar(request.filePath)
        error('postdata:validateRequest:BadPath', 'baseDir and filePath must be text.');
    end
    if ~pd_is_text_scalar(request.progressMode)
        error('postdata:validateRequest:BadProgressMode', ...
            'progressMode must be text.');
    end
    if isempty(strtrim(request.baseDir)) || ~exist(request.baseDir, 'dir')
        error('postdata:validateRequest:BadBaseDirectory', ...
            'baseDir must be an existing directory: %s', request.baseDir);
    end
    request.progressMode = lower(strtrim(request.progressMode));
    if ~any(strcmp(request.progressMode, {'auto','console','waitbar','off'}))
        error('postdata:validateRequest:BadProgressMode', ...
            'progressMode must be auto, console, waitbar, or off.');
    end
    if ~islogical(request.makePlots) || ~isscalar(request.makePlots)
        error('postdata:validateRequest:BadMakePlots', 'makePlots must be a logical scalar.');
    end
    if ~iscell(request.analysisOptions) || mod(numel(request.analysisOptions), 2) ~= 0
        error('postdata:validateRequest:BadOptions', ...
            'analysisOptions must be a name-value cell array.');
    end
    if ~isfield(request, 'plotOptions'), request.plotOptions = struct(); end
    if ~(isempty(request.plotOptions) || ...
            (isstruct(request.plotOptions) && isscalar(request.plotOptions)) || ...
            (iscell(request.plotOptions) && mod(numel(request.plotOptions), 2) == 0))
        error('postdata:validateRequest:BadPlotOptions', ...
            'plotOptions must be a scalar struct or name-value cell array.');
    end
    if ~isstruct(request.chunk) || ~isscalar(request.chunk)
        error('postdata:validateRequest:BadChunkOptions', ...
            'chunk must be a scalar structure.');
    end
    chunkDefaults = struct('dimension', 'auto', 'variable', 'c_rho', ...
        'dV', [], 'coordScale', 1, 'coordRangeX', [], 'coordRangeY', [], ...
        'coordRangeZ', [], ...
        'gradientVariable', '', 'gradientSmoothLevel', 0, ...
        'strainRateVelocityComponent', 'vz', ...
        'strainRateDensityVariable', 'c_rho', 'plotOptions', {{}});
    chunkNames = fieldnames(chunkDefaults);
    for i = 1:numel(chunkNames)
        if ~isfield(request.chunk, chunkNames{i})
            request.chunk.(chunkNames{i}) = chunkDefaults.(chunkNames{i});
        end
    end
    if ~iscell(request.chunk.plotOptions) || ...
            mod(numel(request.chunk.plotOptions), 2) ~= 0
        error('postdata:validateRequest:BadLegacyPlotOptions', ...
            'chunk.plotOptions must be a name-value cell array.');
    end
    if ~isfield(request, 'resultLevel'), request.resultLevel = 'standard'; end
    if ~pd_is_text_scalar(request.resultLevel)
        error('postdata:validateRequest:BadResultLevel', ...
            'resultLevel must be text.');
    end
    request.resultLevel = lower(strtrim(pd_to_char(request.resultLevel)));
    if ~any(strcmpi(request.resultLevel, {'summary','standard','full'}))
        error('postdata:validateRequest:BadResultLevel', ...
            'resultLevel must be summary, standard, or full.');
    end
    if ~isfield(request, 'execution') || ~isstruct(request.execution) || ...
            ~isscalar(request.execution)
        request.execution = struct();
    end
    if ~isfield(request.execution, 'cancelCallback')
        request.execution.cancelCallback = @() false;
    end
    if ~isfield(request.execution, 'progressCallback')
        request.execution.progressCallback = @(fraction, message) [];
    end
    if ~isfield(request.execution, 'logFile')
        request.execution.logFile = '';
    end
    request.execution.logFile = pd_to_char(request.execution.logFile);
    if ~pd_is_text_scalar(request.execution.logFile)
        error('postdata:validateRequest:BadLogFile', ...
            'execution.logFile must be text.');
    end
    if ~isa(request.execution.cancelCallback, 'function_handle') || ...
            ~isa(request.execution.progressCallback, 'function_handle')
        error('postdata:validateRequest:BadExecutionCallback', ...
            'Execution callbacks must be function handles.');
    end

    if ~isstruct(request.selection) || ~isscalar(request.selection)
        error('postdata:validateRequest:BadSelection', ...
            'selection must be a scalar structure.');
    end
    selectionRequired = {'mode', 'value', 'slurmPath', 'slurmModuleIndex'};
    for i = 1:numel(selectionRequired)
        if ~isfield(request.selection, selectionRequired{i})
            error('postdata:validateRequest:BadSelection', ...
                'selection is missing field "%s".', selectionRequired{i});
        end
    end
    request.selection.mode = pd_to_char(request.selection.mode);
    request.selection.slurmPath = pd_to_char(request.selection.slurmPath);
    if ~pd_is_text_scalar(request.selection.mode) || ...
            ~pd_is_text_scalar(request.selection.slurmPath)
        error('postdata:validateRequest:BadSelection', ...
            'selection.mode and selection.slurmPath must be text.');
    end
    switch lower(strtrim(request.selection.mode))
        case 'index'
            request.selection.mode = 'Index';
        case 'timestep'
            request.selection.mode = 'TimeStep';
        case 'time'
            request.selection.mode = 'Time';
        otherwise
            error('postdata:validateRequest:BadSelectionMode', ...
                'selection.mode must be Index, TimeStep, or Time.');
    end
    if ~(isnumeric(request.selection.value) && isscalar(request.selection.value) && ...
            isfinite(request.selection.value))
        error('postdata:validateRequest:BadSelectionValue', ...
            'selection.value must be a finite numeric scalar.');
    end
    if strcmp(request.selection.mode, 'Index') && ...
            (request.selection.value < 1 || ...
             abs(request.selection.value - round(request.selection.value)) > 1e-12)
        error('postdata:validateRequest:BadSelectionIndex', ...
            'Index selection requires a positive integer value.');
    end
    if ~(isnumeric(request.selection.slurmModuleIndex) && ...
            isscalar(request.selection.slurmModuleIndex) && ...
            isfinite(request.selection.slurmModuleIndex) && ...
            request.selection.slurmModuleIndex >= 1 && ...
            abs(request.selection.slurmModuleIndex - ...
                round(request.selection.slurmModuleIndex)) <= 1e-12)
        error('postdata:validateRequest:BadSlurmModuleIndex', ...
            'selection.slurmModuleIndex must be a positive integer.');
    end
    request.selection.slurmModuleIndex = round(request.selection.slurmModuleIndex);
    if strcmpi(request.selection.mode, 'Time') && isempty(request.selection.slurmPath)
        error('postdata:validateRequest:MissingSlurmPath', ...
            'selection.slurmPath is required for Time selection.');
    end
    if strcmpi(request.selection.mode, 'Time') && ...
            ~exist(request.selection.slurmPath, 'file')
        error('postdata:validateRequest:SlurmFileNotFound', ...
            'selection.slurmPath does not exist: %s', request.selection.slurmPath);
    end

    request = normalizeCatalogOptions(request);
end

function request = normalizeCatalogOptions(request)
    catalog = pd_option_catalog(request.analysisType);
    data = pd_catalog_table_data(catalog);
    for i = 1:numel(catalog)
        if strcmp(catalog(i).Target, 'analysis'), continue; end
        data{i, 2} = requestValue(request, catalog(i).Target);
    end

    optionNames = cell(1, numel(request.analysisOptions) / 2);
    for i = 1:2:numel(request.analysisOptions)
        name = request.analysisOptions{i};
        if ~pd_is_text_scalar(name)
            error('postdata:validateRequest:BadOptionName', ...
                'Every analysis option name must be text.');
        end
        name = pd_to_char(name);
        row = find(strcmpi(name, {catalog.Name}), 1, 'first');
        if isempty(row) || ~strcmp(catalog(row).Target, 'analysis')
            error('postdata:validateRequest:UnknownAnalysisOption', ...
                'Unknown calculation option for %s: %s.', ...
                request.analysisType, name);
        end
        optionIndex = (i + 1) / 2;
        optionNames{optionIndex} = lower(catalog(row).Name);
        if nnz(strcmp(optionNames{optionIndex}, optionNames(1:optionIndex))) > 1
            error('postdata:validateRequest:DuplicateAnalysisOption', ...
                'Calculation option is specified more than once: %s.', name);
        end
        data{row, 2} = request.analysisOptions{i + 1};
    end
    request = pd_apply_analysis_options(request, catalog, data);
end

function value = requestValue(request, target)
    switch target
        case 'chunk.dimension'
            value = request.chunk.dimension;
        case 'chunk.variable'
            value = request.chunk.variable;
        case 'chunk.dV'
            value = request.chunk.dV;
        case 'chunk.coordScale'
            value = request.chunk.coordScale;
        case 'chunk.coordRangeX'
            value = request.chunk.coordRangeX;
        case 'chunk.coordRangeY'
            value = request.chunk.coordRangeY;
        case 'chunk.coordRangeZ'
            value = request.chunk.coordRangeZ;
        case 'chunk.gradientVariable'
            value = request.chunk.gradientVariable;
        case 'chunk.gradientSmoothLevel'
            value = request.chunk.gradientSmoothLevel;
        case 'chunk.strainRateVelocityComponent'
            value = request.chunk.strainRateVelocityComponent;
        case 'chunk.strainRateDensityVariable'
            value = request.chunk.strainRateDensityVariable;
        otherwise
            error('postdata:validateRequest:UnknownOptionTarget', ...
                'Unsupported request option target: %s.', target);
    end
end
