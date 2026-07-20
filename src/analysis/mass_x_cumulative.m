function out = mass_x_cumulative(chunkPath, varargin)
%MASS_X_CUMULATIVE SPH mass-x from chunk particle counts.
% Supports planar 2D and volumetric 3D SPH. Chunk Ncount is normalized by
% the transverse y width (and z width for 3D) and accumulated from large x
% to small x.

    p = inputParser;
    p.addRequired('chunkPath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);
    p.addParameter('CancelCallback', @() false, ...
        @(x) isa(x, 'function_handle'));
    p.addParameter('ProgressCallback', @(fraction, message) [], ...
        @(x) isa(x, 'function_handle'));
    p.addParameter('ChunkDim', 'auto', @isTextScalar);
    p.addParameter('CoordinateVar', 'auto', @isTextScalar);
    p.addParameter('CoordinateVarY', 'auto', @isTextScalar);
    p.addParameter('NcountVar', 'Ncount', @isTextScalar);
    p.addParameter('SphDimension', 2, @isnumeric);
    p.addParameter('InitialDensity', [], @isnumeric);
    p.addParameter('ParticleSpacing', [], @isnumeric);
    p.addParameter('RawLengthUnitUm', 10, @isnumeric);
    p.addParameter('UseFileUnitMetadata', true, @islogical);
    p.addParameter('TransverseWidth', [], @isnumeric);
    p.addParameter('OutOfPlaneWidth', [], @isnumeric);
    p.addParameter('SliceCentersY', [], @isnumeric);
    p.addParameter('SliceWidthsY', [], @isnumeric);
    p.addParameter('CoordinateFactor', 10, @isnumeric);
    p.addParameter('CoordinateRange', [], @isnumeric);
    p.addParameter('CoordinateLabel', 'x', @isTextScalar);
    p.addParameter('CoordinateUnit', 'um', @isTextScalar);
    p.addParameter('CumulativeDirection', 'high-to-low', @isTextScalar);
    p.addParameter('MakePlot', true, @islogical);
    p.parse(chunkPath, varargin{:});
    opt = p.Results;
    validateOptions(opt);
    chunkPath = toChar(chunkPath);

    selectorArgs = {'SelectBy', toChar(opt.SelectBy), ...
        'Index', opt.Index, 'TimeStep', opt.TimeStep, 'Time', opt.Time, ...
        'SlurmPath', toChar(opt.SlurmPath), ...
        'SlurmModuleIndex', opt.SlurmModuleIndex, ...
        'ProgressMode', toChar(opt.ProgressMode), ...
        'CancelCallback', opt.CancelCallback, ...
        'ProgressCallback', opt.ProgressCallback};
    step = read_chunk_step_fast(chunkPath, selectorArgs{:});
    units = pd_spid_unit_info(step.unitSystem);
    usedFileUnitMetadata = opt.UseFileUnitMetadata && units.isKnown;
    if usedFileUnitMetadata
        opt.RawLengthUnitUm = units.lengthUmPerUnit;
        opt.CoordinateFactor = metadataCoordinateFactor( ...
            units, opt.CoordinateUnit);
    end

    xVar = resolveVariable(step.colIndex, opt.CoordinateVar, ...
        {'Coord1','c_x','Chunk'}, 'mass_x_cumulative:MissingCoordinateX');
    countVar = resolveVariable(step.colIndex, opt.NcountVar, ...
        {}, 'mass_x_cumulative:MissingNcount');
    dimension = resolveDimension(opt.ChunkDim, opt.CoordinateVarY, step.colIndex);
    xRaw = step.data(:, step.colIndex.(xVar));
    count = step.data(:, step.colIndex.(countVar));
    validateChunkValues(xRaw, count);

    yVar = '';
    if strcmp(dimension, '2d')
        yVar = resolveVariable(step.colIndex, opt.CoordinateVarY, ...
            {'Coord2','c_y'}, 'mass_x_cumulative:MissingCoordinateY');
        yRaw = step.data(:, step.colIndex.(yVar));
        if any(~isfinite(yRaw))
            error('mass_x_cumulative:BadCoordinateY', ...
                'The y coordinate column must contain finite values.');
        end
        [xUniqueRaw, localCount, sliceDefinitions] = ...
            buildTwoDimensionalCounts(xRaw, yRaw, count, opt);
    else
        [xUniqueRaw, localCount] = aggregateByX(xRaw, count);
        sliceDefinitions = oneDimensionalSlice(opt);
    end

    coordinate = xUniqueRaw .* opt.CoordinateFactor;
    if ~isempty(opt.CoordinateRange)
        keep = coordinate >= opt.CoordinateRange(1) & ...
            coordinate <= opt.CoordinateRange(2);
        coordinate = coordinate(keep);
        xUniqueRaw = xUniqueRaw(keep);
        localCount = localCount(keep, :);
    end
    if isempty(coordinate)
        error('mass_x_cumulative:EmptyCoordinateRange', ...
            'CoordinateRange excludes every x bin.');
    end

    rawLengthCm = opt.RawLengthUnitUm .* 1e-4;
    particleMassModelValue = opt.InitialDensity .* ...
        opt.ParticleSpacing .^ opt.SphDimension;
    particleLineMass = NaN;
    particleMass = NaN;
    if opt.SphDimension == 2
        massModel = 'planar-2d-count';
        particleLineMass = opt.InitialDensity .* ...
            (opt.ParticleSpacing .* rawLengthCm) .^ 2;
    else
        massModel = 'volumetric-3d-count';
        particleMass = opt.InitialDensity .* ...
            (opt.ParticleSpacing .* rawLengthCm) .^ 3;
    end
    localDensity = zeros(size(localCount));
    for i = 1:numel(sliceDefinitions)
        widthCm = sliceDefinitions(i).coveredWidthRaw .* rawLengthCm;
        if opt.SphDimension == 2
            localDensity(:, i) = localCount(:, i) .* particleLineMass ./ ...
                widthCm .* 1000;
            sliceDefinitions(i).normalizationAreaCm2 = NaN;
        else
            areaCm2 = widthCm .* opt.OutOfPlaneWidth .* rawLengthCm;
            localDensity(:, i) = localCount(:, i) .* particleMass ./ ...
                areaCm2 .* 1000;
            sliceDefinitions(i).normalizationAreaCm2 = areaCm2;
        end
    end
    cumulativeDensity = flipud(cumsum(flipud(localDensity), 1));
    monotonicTolerance = max(1, max(abs(cumulativeDensity(:)))) .* 1e-12;
    isMonotonic = all(diff(cumulativeDensity, 1, 1) <= ...
        monotonicTolerance, 1);
    labels = {sliceDefinitions.label};
    plotMask = any(localDensity > 0, 1);

    fig = [];
    if opt.MakePlot
        fig = figure('Name', 'SPH Cumulative Areal Density vs Position');
        ax = axes('Parent', fig);
        renderDistribution(ax, coordinate, cumulativeDensity, labels, ...
            opt, step.timestep);
    end

    out = struct();
    out.filePath = chunkPath;
    out.selection = buildSelectionInfo(step, opt);
    out.stepIndex = step.stepIndex;
    out.timestep = step.timestep;
    out.simulationType = 'SPH';
    out.sphDimension = opt.SphDimension;
    out.massModel = massModel;
    out.chunkDimension = dimension;
    out.coordinate = coordinate;
    out.x = coordinate;
    out.coordinateRaw = xUniqueRaw;
    out.coordinateVar = xVar;
    out.coordinateVarY = yVar;
    out.countVar = countVar;
    out.coordinateFactor = opt.CoordinateFactor;
    out.coordinateRange = opt.CoordinateRange;
    out.coordinateLabel = toChar(opt.CoordinateLabel);
    out.coordinateUnit = toChar(opt.CoordinateUnit);
    out.initialDensityGcm3 = opt.InitialDensity;
    out.particleSpacingRaw = opt.ParticleSpacing;
    out.rawLengthUnitUm = opt.RawLengthUnitUm;
    out.useFileUnitMetadata = opt.UseFileUnitMetadata;
    out.usedFileUnitMetadata = usedFileUnitMetadata;
    out.unitSystem = step.unitSystem;
    out.physicalTime = step.physicalTime;
    out.inputFormat = step.inputFormat;
    out.taskName = step.taskName;
    out.chunkKind = step.chunkKind;
    out.particleMassModelValue = particleMassModelValue;
    out.particleMassG = particleMass;
    out.particleLineMassGPerCm = particleLineMass;
    out.transverseWidthRaw = opt.TransverseWidth;
    out.outOfPlaneWidthRaw = opt.OutOfPlaneWidth;
    out.sliceCentersYRaw = opt.SliceCentersY(:).';
    out.sliceWidthsYRaw = opt.SliceWidthsY(:).';
    out.sliceBoundaryMethod = 'fractional-bin-overlap';
    out.sliceDefinitions = sliceDefinitions;
    out.localParticleCount = localCount;
    out.density = localDensity;
    out.cumulativeDensity = cumulativeDensity;
    out.densityVars = labels;
    out.plottedVars = labels(plotMask);
    out.cumulativeDirection = 'high-to-low';
    out.isMonotonicDecreasing = isMonotonic;
    out.figure = fig;
end

function [xUnique, localCount, definitions] = ...
        buildTwoDimensionalCounts(xRaw, yRaw, count, opt)

    [xUnique, ~, xGroup] = unique(xRaw);
    [yUnique, ~, yGroup] = unique(yRaw);
    if numel(yUnique) > 1
        spacing = diff(yUnique);
        dy = median(spacing);
        tolerance = max(1, abs(dy)) .* 1e-8;
        if dy <= 0 || any(abs(spacing - dy) > tolerance)
            error('mass_x_cumulative:NonUniformYGrid', ...
                'SPH mass-x currently requires a uniform Coord2/c_y grid.');
        end
    else
        if isempty(opt.TransverseWidth)
            error('mass_x_cumulative:NeedSingleRowWidth', ...
                ['A 2D chunk with one y row requires TransverseWidth to ', ...
                 'define that row width in raw coordinate units.']);
        end
        dy = opt.TransverseWidth;
    end

    cellLower = yUnique - dy ./ 2;
    cellUpper = yUnique + dy ./ 2;
    domainLower = cellLower(1);
    domainUpper = cellUpper(end);
    [centers, widths, labels, isFull] = resolveSlices( ...
        opt.SliceCentersY, opt.SliceWidthsY, domainLower, domainUpper, ...
        opt.CoordinateFactor, opt.CoordinateUnit);

    localCount = zeros(numel(xUnique), numel(centers));
    definitions = repmat(emptySliceDefinition(), 1, numel(centers));
    for i = 1:numel(centers)
        requestedLower = centers(i) - widths(i) ./ 2;
        requestedUpper = centers(i) + widths(i) ./ 2;
        overlap = max(0, min(cellUpper, requestedUpper) - ...
            max(cellLower, requestedLower));
        coveredWidth = sum(overlap);
        if coveredWidth <= max(1, widths(i)) * eps
            error('mass_x_cumulative:SliceOutsideDomain', ...
                'The y slice centered at %.8g does not overlap the chunk grid.', ...
                centers(i));
        end
        fractions = overlap ./ dy;
        weightedCount = count .* fractions(yGroup);
        localCount(:, i) = accumarray(xGroup, weightedCount, ...
            [numel(xUnique), 1], @sum, 0);

        definitions(i).label = labels{i};
        definitions(i).isFullWidth = isFull(i);
        definitions(i).centerRaw = centers(i);
        definitions(i).centerDisplay = centers(i) .* opt.CoordinateFactor;
        definitions(i).requestedWidthRaw = widths(i);
        definitions(i).requestedWidthDisplay = widths(i) .* opt.CoordinateFactor;
        definitions(i).coveredWidthRaw = coveredWidth;
        definitions(i).coveredWidthDisplay = coveredWidth .* opt.CoordinateFactor;
        definitions(i).lowerRaw = max(requestedLower, domainLower);
        definitions(i).upperRaw = min(requestedUpper, domainUpper);
        definitions(i).gridSpacingYRaw = dy;
    end
end

function [centers, widths, labels, isFull] = resolveSlices( ...
        centerInput, widthInput, domainLower, domainUpper, ...
        coordinateFactor, coordinateUnit)

    domainWidth = domainUpper - domainLower;
    if isempty(centerInput)
        if ~isempty(widthInput)
            error('mass_x_cumulative:SliceWidthWithoutCenter', ...
                'SliceWidthsY requires SliceCentersY.');
        end
        centers = (domainLower + domainUpper) ./ 2;
        widths = domainWidth;
        labels = {'full_y'};
        isFull = true;
        return;
    end

    centers = centerInput(:).';
    widths = widthInput(:).';
    if isempty(widths)
        error('mass_x_cumulative:MissingSliceWidths', ...
            'SliceWidthsY is required when SliceCentersY is specified.');
    end
    if isscalar(widths)
        widths = repmat(widths, size(centers));
    elseif numel(widths) ~= numel(centers)
        error('mass_x_cumulative:SliceSizeMismatch', ...
            'SliceWidthsY must be scalar or match SliceCentersY in length.');
    end
    if any(~isfinite(centers)) || any(~isfinite(widths)) || any(widths <= 0)
        error('mass_x_cumulative:BadSlices', ...
            'Slice centers must be finite and slice widths must be positive.');
    end

    labels = cell(1, numel(centers));
    isFull = false(size(centers));
    unit = regexprep(strtrim(toChar(coordinateUnit)), '[^A-Za-z0-9]', '');
    if ~isempty(unit)
        unit = ['_', unit];
    end
    for i = 1:numel(centers)
        labels{i} = sprintf('y=%s%s_width=%s%s', ...
            formatNumber(centers(i) .* coordinateFactor), unit, ...
            formatNumber(widths(i) .* coordinateFactor), unit);
    end
end

function definition = oneDimensionalSlice(opt)
    width = opt.TransverseWidth;
    if isempty(width)
        error('mass_x_cumulative:NeedTransverseWidth', ...
            ['A bin1d file contains no y extent. Specify TransverseWidth ', ...
             'in raw SPH coordinate units.']);
    end
    definition = emptySliceDefinition();
    definition.label = 'full_y';
    definition.isFullWidth = true;
    definition.centerRaw = NaN;
    definition.centerDisplay = NaN;
    definition.requestedWidthRaw = width;
    definition.requestedWidthDisplay = width .* opt.CoordinateFactor;
    definition.coveredWidthRaw = width;
    definition.coveredWidthDisplay = width .* opt.CoordinateFactor;
    definition.lowerRaw = NaN;
    definition.upperRaw = NaN;
    definition.gridSpacingYRaw = NaN;
end

function definition = emptySliceDefinition()
    definition = struct('label', '', 'isFullWidth', false, ...
        'centerRaw', NaN, 'centerDisplay', NaN, ...
        'requestedWidthRaw', NaN, 'requestedWidthDisplay', NaN, ...
        'coveredWidthRaw', NaN, 'coveredWidthDisplay', NaN, ...
        'lowerRaw', NaN, 'upperRaw', NaN, ...
        'gridSpacingYRaw', NaN, 'normalizationAreaCm2', NaN);
end

function [xUnique, total] = aggregateByX(xRaw, values)
    [xUnique, ~, groups] = unique(xRaw);
    total = accumarray(groups, values, [numel(xUnique), 1], @sum, 0);
end

function dimension = resolveDimension(requested, requestedY, col)
    requested = lower(strtrim(toChar(requested)));
    hasY = isfield(col, 'Coord2') || isfield(col, 'c_y');
    requestedY = strtrim(toChar(requestedY));
    if ~strcmpi(requestedY, 'auto')
        hasY = hasY || isfield(col, matlab.lang.makeValidName(requestedY));
    end
    switch requested
        case 'auto'
            if hasY
                dimension = '2d';
            else
                dimension = '1d';
            end
        case '1d'
            dimension = '1d';
        case '2d'
            if ~hasY
                error('mass_x_cumulative:MissingCoordinateY', ...
                    'ChunkDim=2d requires Coord2 or c_y.');
            end
            dimension = '2d';
        otherwise
            error('mass_x_cumulative:BadChunkDim', ...
                'ChunkDim must be auto, 1d, or 2d.');
    end
end

function name = resolveVariable(col, requested, autoCandidates, errorId)
    requested = strtrim(toChar(requested));
    if strcmpi(requested, 'auto')
        candidates = autoCandidates;
    else
        candidates = {matlab.lang.makeValidName(requested)};
    end
    name = '';
    for i = 1:numel(candidates)
        if isfield(col, candidates{i})
            name = candidates{i};
            break;
        end
    end
    if isempty(name)
        if isempty(candidates)
            expected = requested;
        else
            expected = strjoin(candidates, ', ');
        end
        error(errorId, 'Required chunk column not found. Expected: %s.', expected);
    end
end

function validateChunkValues(xRaw, count)
    if any(~isfinite(xRaw))
        error('mass_x_cumulative:BadCoordinateX', ...
            'The x coordinate column must contain finite values.');
    end
    if any(~isfinite(count)) || any(count < 0)
        error('mass_x_cumulative:BadNcount', ...
            'Ncount must contain finite non-negative particle counts.');
    end
end

function renderDistribution(ax, coordinate, cumulativeDensity, labels, opt, timestep)
    if isempty(cumulativeDensity) || ~any(cumulativeDensity(:) > 0)
        text(0.5, 0.5, 'No particles in selected region', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
    else
        handles = plot(ax, coordinate, cumulativeDensity, 'LineWidth', 1.6);
        for i = 1:numel(handles)
            set(handles(i), 'DisplayName', labels{i});
        end
        legend(ax, 'show', 'Location', 'best', 'Interpreter', 'none');
    end
    xlabel(ax, coordinateAxisLabel(opt));
    ylabel(ax, 'Cumulative Areal Density (mg/cm^2)');
    title(ax, sprintf('SPH mass-x cumulative @ %g', timestep));
    grid(ax, 'on');
end

function label = coordinateAxisLabel(opt)
    label = toChar(opt.CoordinateLabel);
    unit = strtrim(toChar(opt.CoordinateUnit));
    if ~isempty(unit)
        label = sprintf('%s (%s)', label, unit);
    end
end

function validateOptions(opt)
    validatePositiveScalar(opt.InitialDensity, 'InitialDensity');
    validatePositiveScalar(opt.ParticleSpacing, 'ParticleSpacing');
    validatePositiveScalar(opt.RawLengthUnitUm, 'RawLengthUnitUm');
    validatePositiveScalar(opt.CoordinateFactor, 'CoordinateFactor');
    if ~(isscalar(opt.SphDimension) && isfinite(opt.SphDimension) && ...
            any(opt.SphDimension == [2, 3]))
        error('mass_x_cumulative:BadSphDimension', ...
            'SphDimension must be 2 or 3.');
    end
    if ~isempty(opt.TransverseWidth)
        validatePositiveScalar(opt.TransverseWidth, 'TransverseWidth');
    end
    if opt.SphDimension == 3
        validatePositiveScalar(opt.OutOfPlaneWidth, 'OutOfPlaneWidth');
    elseif ~isempty(opt.OutOfPlaneWidth)
        validatePositiveScalar(opt.OutOfPlaneWidth, 'OutOfPlaneWidth');
    end
    if ~isempty(opt.CoordinateRange) && ~(isnumeric(opt.CoordinateRange) && ...
            numel(opt.CoordinateRange) == 2 && ...
            all(isfinite(opt.CoordinateRange)) && ...
            opt.CoordinateRange(2) >= opt.CoordinateRange(1))
        error('mass_x_cumulative:BadCoordinateRange', ...
            'CoordinateRange must be empty or [min max].');
    end
    if ~strcmpi(strtrim(toChar(opt.CumulativeDirection)), 'high-to-low')
        error('mass_x_cumulative:DirectionMustDecrease', ...
            ['mass-x is defined as accumulation from large x to small x; ', ...
             'CumulativeDirection must be high-to-low.']);
    end
end

function factor = metadataCoordinateFactor(units, requestedUnit)
    unit = lower(regexprep(strtrim(toChar(requestedUnit)), '\s+', ''));
    switch unit
        case {'um','micrometer','micrometers'}
            factor = units.lengthUmPerUnit;
        case {'nm','nanometer','nanometers'}
            factor = 1000 .* units.lengthUmPerUnit;
        case {'mm','millimeter','millimeters'}
            factor = 1e-3 .* units.lengthUmPerUnit;
        case {'cm','centimeter','centimeters'}
            factor = 1e-4 .* units.lengthUmPerUnit;
        case {'m','meter','meters'}
            factor = 1e-6 .* units.lengthUmPerUnit;
        case {'native','raw'}
            factor = 1;
        otherwise
            error('mass_x_cumulative:UnsupportedCoordinateUnit', ...
                ['CoordinateUnit "%s" cannot be derived from SPID Units ' ...
                 'metadata. Use nm, um, mm, cm, m, native, or set ' ...
                 'UseFileUnitMetadata=false with a manual CoordinateFactor.'], ...
                toChar(requestedUnit));
    end
end

function validatePositiveScalar(value, name)
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
        error('mass_x_cumulative:BadPhysicalOption', ...
            '%s must be a positive finite scalar.', name);
    end
end

function info = buildSelectionInfo(step, opt)
    mode = lower(strtrim(toChar(opt.SelectBy)));
    info = struct('mode', mode, 'stepIndex', step.stepIndex, ...
        'mappedTimeStep', step.timestep);
    switch mode
        case 'index'
            info.requestedIndex = opt.Index;
        case 'timestep'
            info.requestedTimeStep = opt.TimeStep;
        case 'time'
            info.requestedTime = opt.Time;
    end
end

function text = formatNumber(value)
    text = sprintf('%.8g', value);
end

function tf = isTextScalar(v)
    tf = ischar(v) || (isstring(v) && isscalar(v));
end

function s = toChar(v)
    if isstring(v)
        s = char(v);
    else
        s = v;
    end
end
