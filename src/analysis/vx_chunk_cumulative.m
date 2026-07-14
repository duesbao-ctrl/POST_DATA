function out = vx_chunk_cumulative(vxChunkPath, varargin)
%VX_CHUNK_CUMULATIVE Cumulative areal density versus velocity.

    p = inputParser;
    p.addRequired('vxChunkPath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);
    p.addParameter('CancelCallback', @() false, @(x) isa(x, 'function_handle'));
    p.addParameter('ProgressCallback', @(fraction, message) [], ...
        @(x) isa(x, 'function_handle'));
    p.addParameter('VelocityVar', 'Chunk', @isTextScalar);
    p.addParameter('VelocityFactor', 0.001, @isnumeric);
    p.addParameter('VelocityLabel', 'Velocity', @isTextScalar);
    p.addParameter('VelocityUnit', 'km/s', @isTextScalar);
    p.addParameter('CumulativeDirection', 'high-to-low', @isTextScalar);
    p.addParameter('NegativeDensityMode', 'clip', @isTextScalar);
    p.addParameter('DensityVars', {}, @(x) iscell(x) || isTextScalar(x));
    p.addParameter('MakePlot', true, @islogical);
    p.parse(vxChunkPath, varargin{:});
    opt = p.Results;
    validatePositiveFactor(opt.VelocityFactor);
    vxChunkPath = toChar(vxChunkPath);

    selectorArgs = {'SelectBy', toChar(opt.SelectBy), 'Index', opt.Index, ...
        'TimeStep', opt.TimeStep, 'Time', opt.Time, ...
        'SlurmPath', toChar(opt.SlurmPath), ...
        'SlurmModuleIndex', opt.SlurmModuleIndex, ...
        'ProgressMode', toChar(opt.ProgressMode), ...
        'CancelCallback', opt.CancelCallback, ...
        'ProgressCallback', opt.ProgressCallback};
    step = read_chunk_step_fast(vxChunkPath, selectorArgs{:});
    velocityVar = matlab.lang.makeValidName(toChar(opt.VelocityVar));
    if ~isfield(step.colIndex, velocityVar)
        error('vx_chunk_cumulative:MissingVelocityVar', ...
            'Velocity coordinate column "%s" is required.', velocityVar);
    end
    velocity = step.data(:, step.colIndex.(velocityVar)) .* opt.VelocityFactor;
    distribution = pd_cumulative_areal_density(step, velocity, ...
        opt.DensityVars, opt.CumulativeDirection, 'vx_chunk_cumulative', ...
        opt.NegativeDensityMode);

    fig = [];
    if opt.MakePlot
        fig = figure('Name', 'Cumulative Areal Density vs Velocity');
        ax = axes('Parent', fig);
        renderDistribution(ax, distribution, step.timestep, opt);
    end

    out = struct();
    out.filePath = vxChunkPath;
    out.selection = buildSelectionInfo(step, opt);
    out.stepIndex = step.stepIndex;
    out.timestep = step.timestep;
    out.velocity = distribution.coordinate;
    out.density = distribution.density;
    out.cumulativeDensity = distribution.cumulativeDensity;
    out.densityVars = distribution.densityVars;
    out.plottedVars = distribution.plottedVars;
    out.cumulativeDirection = distribution.direction;
    out.negativeDensityMode = distribution.negativeDensityMode;
    out.negativeValueCount = distribution.negativeValueCount;
    out.isMonotonicDecreasing = distribution.isMonotonicDecreasing;
    out.velocityFactor = opt.VelocityFactor;
    out.velocityVar = velocityVar;
    out.velocityLabel = toChar(opt.VelocityLabel);
    out.velocityUnit = toChar(opt.VelocityUnit);
    out.figure = fig;
end

function renderDistribution(ax, distribution, timestep, opt)
    indices = find(distribution.plotMask);
    if isempty(indices)
        text(0.5, 0.5, 'No non-zero cumulative curves', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
    else
        handles = plot(ax, distribution.coordinate, ...
            distribution.cumulativeDensity(:, indices), 'LineWidth', 1.6);
        for i = 1:numel(handles)
            set(handles(i), 'DisplayName', distribution.densityVars{indices(i)});
        end
        legend(ax, 'show', 'Location', 'best', 'Interpreter', 'none');
    end
    xlabel(ax, axisLabel(opt.VelocityLabel, opt.VelocityUnit));
    ylabel(ax, 'Cumulative Areal Density (mg/cm^2)');
    title(ax, sprintf('mass-v cumulative distribution @ timestep %g', timestep));
    grid(ax, 'on');
end

function label = axisLabel(name, unit)
    label = toChar(name);
    unit = strtrim(toChar(unit));
    if ~isempty(unit), label = sprintf('%s (%s)', label, unit); end
end

function info = buildSelectionInfo(step, opt)
    mode = lower(strtrim(toChar(opt.SelectBy)));
    info = struct('mode', mode, 'stepIndex', step.stepIndex, ...
        'mappedTimeStep', step.timestep);
    switch mode
        case 'index', info.requestedIndex = opt.Index;
        case 'timestep', info.requestedTimeStep = opt.TimeStep;
        case 'time', info.requestedTime = opt.Time;
    end
end

function validatePositiveFactor(value)
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
        error('vx_chunk_cumulative:BadVelocityFactor', ...
            'VelocityFactor must be a positive finite scalar.');
    end
end

function tf = isTextScalar(v)
    tf = ischar(v) || (isstring(v) && isscalar(v));
end

function s = toChar(v)
    if isstring(v), s = char(v); else, s = v; end
end
