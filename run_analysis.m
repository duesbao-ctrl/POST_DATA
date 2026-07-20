function result = run_analysis(taskType, varargin)
%RUN_ANALYSIS Compatibility wrapper around the modular POST_DATA request API.
%   New scripts should prefer pd_create_request + postdata_run. This wrapper
%   preserves the established chunk/cluster/vx/network2d calling convention
%   and adds mass-v and mass-x aliases without duplicating analysis logic.

    postdata_startup();
    p = inputParser;
    p.addRequired('taskType', @pd_is_text_scalar);
    p.addParameter('BaseDir', pwd, @pd_is_text_scalar);
    p.addParameter('SelectBy', 'Index', @pd_is_text_scalar);
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @pd_is_text_scalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @pd_is_text_scalar);
    p.addParameter('ResultLevel', 'full', @pd_is_text_scalar);
    p.addParameter('LogFile', '', @pd_is_text_scalar);
    p.addParameter('PublicationPlotOptions', struct(), @isPlotOptions);

    p.addParameter('ChunkDim', 'auto', @pd_is_text_scalar);
    p.addParameter('ChunkFile', '', @pd_is_text_scalar);
    p.addParameter('Variable', 'c_rho', @pd_is_text_scalar);
    p.addParameter('dV', [], @isnumeric);
    p.addParameter('DoPlot', true, @islogical);
    p.addParameter('PlotOptions', {}, @iscell);
    p.addParameter('CoordScale', 1, @isnumeric);
    p.addParameter('CoordRangeX', [], @isnumeric);
    p.addParameter('CoordRangeY', [], @isnumeric);
    p.addParameter('GradientVariable', '', @pd_is_text_scalar);
    p.addParameter('GradientSmoothLevel', 0, @isnumeric);
    p.addParameter('StrainRateVelocityComponent', 'vz', @pd_is_text_scalar);
    p.addParameter('StrainRateDensityVariable', 'c_rho', @pd_is_text_scalar);

    p.addParameter('ClusterFile', '', @pd_is_text_scalar);
    p.addParameter('ClusterOptions', {}, @iscell);
    p.addParameter('VxFile', '', @pd_is_text_scalar);
    p.addParameter('VxOptions', {}, @iscell);
    p.addParameter('MassXFile', '', @pd_is_text_scalar);
    p.addParameter('MassXOptions', {}, @iscell);
    p.addParameter('NetworkOptions', {}, @iscell);
    p.parse(taskType, varargin{:});
    opt = p.Results;

    type = pd_normalize_analysis_type(taskType);
    request = pd_create_request(type);
    request.baseDir = pd_to_char(opt.BaseDir);
    request.progressMode = pd_to_char(opt.ProgressMode);
    request.resultLevel = pd_to_char(opt.ResultLevel);
    request.plotOptions = opt.PublicationPlotOptions;
    request.execution.logFile = pd_to_char(opt.LogFile);
    request.selection.mode = pd_to_char(opt.SelectBy);
    request.selection.value = selectionValue(opt);
    request.selection.slurmPath = pd_to_char(opt.SlurmPath);
    request.selection.slurmModuleIndex = opt.SlurmModuleIndex;
    switch type
        case 'chunk'
            request.filePath = legacyFile(opt.ChunkFile, 'ChunkFile');
            request.makePlots = opt.DoPlot;
            request.chunk.dimension = pd_to_char(opt.ChunkDim);
            request.chunk.variable = pd_to_char(opt.Variable);
            request.chunk.dV = opt.dV;
            request.chunk.coordScale = opt.CoordScale;
            request.chunk.coordRangeX = opt.CoordRangeX;
            request.chunk.coordRangeY = opt.CoordRangeY;
            request.chunk.gradientVariable = pd_to_char(opt.GradientVariable);
            request.chunk.gradientSmoothLevel = opt.GradientSmoothLevel;
            request.chunk.strainRateVelocityComponent = ...
                pd_to_char(opt.StrainRateVelocityComponent);
            request.chunk.strainRateDensityVariable = ...
                pd_to_char(opt.StrainRateDensityVariable);
            request.chunk.plotOptions = opt.PlotOptions;
        case 'cluster'
            request.filePath = legacyFile(opt.ClusterFile, 'ClusterFile');
            [request.analysisOptions, request.makePlots] = ...
                extractLogicalOption(opt.ClusterOptions, 'MakePlots', true);
        case 'vx'
            request.filePath = legacyFile(opt.VxFile, 'VxFile');
            [request.analysisOptions, request.makePlots] = ...
                extractLogicalOption(opt.VxOptions, 'MakePlot', true);
        case 'massx'
            request.filePath = legacyFile(opt.MassXFile, 'MassXFile');
            [request.analysisOptions, request.makePlots] = ...
                extractLogicalOption(opt.MassXOptions, 'MakePlot', true);
        case 'network2d'
            if strcmpi(pd_to_char(opt.ChunkDim), '1d')
                error('run_analysis:BadNetworkChunkDim', ...
                    'taskType=''network2d'' requires 2D input.');
            end
            request.filePath = legacyFile(opt.ChunkFile, 'ChunkFile');
            request.chunk.coordScale = opt.CoordScale;
            [request.analysisOptions, request.makePlots] = ...
                extractLogicalOption(opt.NetworkOptions, 'MakePlots', true);
    end
    result = postdata_run(request);
end

function value = selectionValue(opt)
    switch lower(strtrim(pd_to_char(opt.SelectBy)))
        case 'index', value = opt.Index;
        case 'timestep', value = opt.TimeStep;
        case 'time', value = opt.Time;
        otherwise
            error('run_analysis:BadSelectBy', ...
                'SelectBy must be Index, TimeStep, or Time.');
    end
    if ~(isnumeric(value) && isscalar(value) && isfinite(value))
        error('run_analysis:MissingSelectionValue', ...
            'The selected Index, TimeStep, or Time value must be a finite scalar.');
    end
end

function fileName = legacyFile(value, argumentName)
    fileName = pd_to_char(value);
    if isempty(fileName), return; end
    if isempty(regexpi(fileName, '\.txt$', 'once'))
        error('run_analysis:MissingTxtSuffix', ...
            '%s must include the .txt suffix to avoid matching cache files: %s', ...
            argumentName, fileName);
    end
end

function [remaining, value] = extractLogicalOption(options, name, defaultValue)
    if mod(numel(options), 2) ~= 0
        error('run_analysis:BadOptions', ...
            'Legacy calculation options must be name-value pairs.');
    end
    value = defaultValue;
    remaining = {};
    found = false;
    for i = 1:2:numel(options)
        if ~pd_is_text_scalar(options{i})
            error('run_analysis:BadOptionName', ...
                'Legacy calculation option names must be text.');
        end
        if strcmpi(pd_to_char(options{i}), name)
            if found
                error('run_analysis:DuplicatePlotFlag', ...
                    '%s is specified more than once.', name);
            end
            value = options{i + 1};
            found = true;
        else
            remaining(end + 1:end + 2) = options(i:i + 1); %#ok<AGROW>
        end
    end
    if ~(islogical(value) && isscalar(value))
        error('run_analysis:BadPlotFlag', '%s must be a logical scalar.', name);
    end
end

function tf = isPlotOptions(value)
    tf = isempty(value) || (isstruct(value) && isscalar(value)) || ...
        (iscell(value) && mod(numel(value), 2) == 0);
end
