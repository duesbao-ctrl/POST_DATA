function out = mass_x_cumulative(chunkPath, varargin)
%MASS_X_CUMULATIVE Cumulative areal density versus spatial coordinate.
% Default accumulation is from the largest coordinate to the smallest.

    p = inputParser;
    p.addRequired('chunkPath', @isTextScalar);
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
    p.addParameter('CoordinateVar', 'auto', @isTextScalar);
    p.addParameter('CoordinateFactor', 1, @isnumeric);
    p.addParameter('CoordinateRange', [], @isnumeric);
    p.addParameter('CoordinateLabel', 'x', @isTextScalar);
    p.addParameter('CoordinateUnit', 'cm', @isTextScalar);
    p.addParameter('CumulativeDirection', 'high-to-low', @isTextScalar);
    p.addParameter('NegativeDensityMode', 'clip', @isTextScalar);
    p.addParameter('DensityVars', {}, @(x) iscell(x) || isTextScalar(x));
    p.addParameter('MakePlot', true, @islogical);
    p.parse(chunkPath, varargin{:});
    opt = p.Results;
    validateOptions(opt);
    chunkPath = toChar(chunkPath);

    selectorArgs = {'SelectBy', toChar(opt.SelectBy), 'Index', opt.Index, ...
        'TimeStep', opt.TimeStep, 'Time', opt.Time, ...
        'SlurmPath', toChar(opt.SlurmPath), ...
        'SlurmModuleIndex', opt.SlurmModuleIndex, ...
        'ProgressMode', toChar(opt.ProgressMode), ...
        'CancelCallback', opt.CancelCallback, ...
        'ProgressCallback', opt.ProgressCallback};
    step = read_chunk_step_fast(chunkPath, selectorArgs{:});
    coordinateVar = resolveCoordinateVariable(step.colIndex, opt.CoordinateVar);
    coordinate = step.data(:, step.colIndex.(coordinateVar)) .* opt.CoordinateFactor;
    if ~isempty(opt.CoordinateRange)
        inRange = coordinate >= opt.CoordinateRange(1) & ...
            coordinate <= opt.CoordinateRange(2);
        step.data = step.data(inRange, :);
        coordinate = coordinate(inRange);
    end
    distribution = pd_cumulative_areal_density(step, coordinate, ...
        opt.DensityVars, opt.CumulativeDirection, 'mass_x_cumulative', ...
        opt.NegativeDensityMode);

    fig = [];
    if opt.MakePlot
        fig = figure('Name', 'Cumulative Areal Density vs Position');
        ax = axes('Parent', fig);
        renderDistribution(ax, distribution, opt, step.timestep);
    end

    out = struct();
    out.filePath = chunkPath;
    out.selection = buildSelectionInfo(step, opt);
    out.stepIndex = step.stepIndex;
    out.timestep = step.timestep;
    out.coordinate = distribution.coordinate;
    out.x = distribution.coordinate;
    out.coordinateVar = coordinateVar;
    out.coordinateFactor = opt.CoordinateFactor;
    out.coordinateRange = opt.CoordinateRange;
    out.coordinateLabel = toChar(opt.CoordinateLabel);
    out.coordinateUnit = toChar(opt.CoordinateUnit);
    out.density = distribution.density;
    out.cumulativeDensity = distribution.cumulativeDensity;
    out.densityVars = distribution.densityVars;
    out.plottedVars = distribution.plottedVars;
    out.cumulativeDirection = distribution.direction;
    out.negativeDensityMode = distribution.negativeDensityMode;
    out.negativeValueCount = distribution.negativeValueCount;
    out.isMonotonicDecreasing = distribution.isMonotonicDecreasing;
    out.figure = fig;
end

function name = resolveCoordinateVariable(col, requested)
    requested = strtrim(toChar(requested));
    if strcmpi(requested, 'auto')
        candidates = {'Coord1','c_x','Chunk'};
    else
        candidates = {matlab.lang.makeValidName(requested)};
    end
    name = '';
    for i = 1:numel(candidates)
        if isfield(col, candidates{i}), name = candidates{i}; break; end
    end
    if isempty(name)
        error('mass_x_cumulative:MissingCoordinate', ...
            'Coordinate column not found. Requested: %s.', requested);
    end
end

function renderDistribution(ax, distribution, opt, timestep)
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
    xlabel(ax, coordinateAxisLabel(opt));
    ylabel(ax, 'Cumulative Areal Density (mg/cm^2)');
    title(ax, sprintf('mass-x cumulative distribution @ timestep %g', timestep));
    grid(ax, 'on');
end

function label = coordinateAxisLabel(opt)
    label = toChar(opt.CoordinateLabel);
    unit = strtrim(toChar(opt.CoordinateUnit));
    if ~isempty(unit), label = sprintf('%s (%s)', label, unit); end
end

function validateOptions(opt)
    if ~(isscalar(opt.CoordinateFactor) && isfinite(opt.CoordinateFactor) && opt.CoordinateFactor > 0)
        error('mass_x_cumulative:BadCoordinateFactor', ...
            'CoordinateFactor must be a positive finite scalar.');
    end
    if ~isempty(opt.CoordinateRange) && ~(isnumeric(opt.CoordinateRange) && ...
            numel(opt.CoordinateRange) == 2 && all(isfinite(opt.CoordinateRange)) && ...
            opt.CoordinateRange(2) >= opt.CoordinateRange(1))
        error('mass_x_cumulative:BadCoordinateRange', ...
            'CoordinateRange must be empty or [min max].');
    end
    if ~strcmpi(strtrim(toChar(opt.CumulativeDirection)), 'high-to-low')
        error('mass_x_cumulative:DirectionMustDecrease', ...
            ['mass-x is defined as accumulation from large coordinate to ', ...
             'small coordinate; CumulativeDirection must be high-to-low.']);
    end
    mode = lower(strtrim(toChar(opt.NegativeDensityMode)));
    if ~any(strcmp(mode, {'clip','error'}))
        error('mass_x_cumulative:BadNegativeDensityMode', ...
            'NegativeDensityMode must be clip or error.');
    end
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

function tf = isTextScalar(v)
    tf = ischar(v) || (isstring(v) && isscalar(v));
end

function s = toChar(v)
    if isstring(v), s = char(v); else, s = v; end
end
