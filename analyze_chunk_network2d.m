function out = analyze_chunk_network2d(chunkFile, varargin)
%ANALYZE_CHUNK_NETWORK2D Analyze pore/matrix network structure on 2D chunk data.
%   out = ANALYZE_CHUNK_NETWORK2D(chunkFile, ...)
%
% A cell is treated as pore when Ncount < ThresholdN, otherwise matrix.
% Connected components are identified on the 2D grid with configurable
% boundary conditions. The function returns global morphology metrics,
% per-phase topology/statistics, and optional plots.

    p = inputParser;
    p.addRequired('chunkFile', @isTextScalar);

    p.addParameter('SelectBy', 'Index', @isTextScalar);
    p.addParameter('Index', 1, @isnumeric);
    p.addParameter('TimeStep', [], @isnumeric);
    p.addParameter('Time', [], @isnumeric);
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);

    p.addParameter('ThresholdN', [], @isnumeric);
    p.addParameter('Connectivity', 4, @isnumeric);
    p.addParameter('Boundary', 'open', @isTextScalar);
    p.addParameter('Dx', [], @isnumeric);
    p.addParameter('Dy', [], @isnumeric);
    p.addParameter('CoordScale', 1, @isnumeric);
    p.addParameter('NcountVar', 'Ncount', @isTextScalar);
    p.addParameter('PositionAxis', 'x', @isTextScalar);
    p.addParameter('PositionNumBins', 20, @isnumeric);
    p.addParameter('PositionRangeX', [], @isnumeric);
    p.addParameter('PositionRangeY', [], @isnumeric);
    p.addParameter('PlotRangeX', [], @isnumeric);
    p.addParameter('PlotRangeY', [], @isnumeric);
    p.addParameter('ProfileAxis', 'x', @isTextScalar);
    p.addParameter('ProfileNumBins', 20, @isnumeric);
    p.addParameter('ProfileRangeX', [], @isnumeric);
    p.addParameter('ProfileRangeY', [], @isnumeric);
    p.addParameter('HistNumBins', 30, @isnumeric);
    p.addParameter('MeanPowerM', 1, @isnumeric);
    p.addParameter('MeanPowerN', 0, @isnumeric);
    p.addParameter('EnableSkeletonGraph', false, @islogical);
    p.addParameter('EnableEvolution', false, @islogical);
    p.addParameter('EvolutionSelectBy', '', @isTextScalar);
    p.addParameter('EvolutionRange', [], @isnumeric);
    p.addParameter('EvolutionStride', 1, @isnumeric);
    p.addParameter('MakePlots', true, @islogical);
    p.parse(chunkFile, varargin{:});
    opt = p.Results;

    chunkFile = toChar(chunkFile);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);
    opt.Boundary = lower(strtrim(toChar(opt.Boundary)));
    opt.NcountVar = toChar(opt.NcountVar);
    opt.PositionAxis = lower(strtrim(toChar(opt.PositionAxis)));
    opt.ProfileAxis = lower(strtrim(toChar(opt.ProfileAxis)));
    opt.EvolutionSelectBy = lower(strtrim(toChar(opt.EvolutionSelectBy)));

    validateInputs(opt);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex, ...
                    'ProgressMode', opt.ProgressMode};

    step = read_chunk_step_fast(chunkFile, selectorArgs{:});
    selectionInfo = buildSelectionInfo(step, opt);
    data = step.data;
    col = step.colIndex;

    [x, y, ncount] = extractRequiredColumns(data, col, opt.NcountVar);
    [x, y] = scaleCoordinates(x, y, opt.CoordScale);
    [ncountGrid, validMask, xCenters, yCenters, presentCount] = buildGrid(x, y, ncount);
    [dx, dy, spacingSource] = resolveGridSpacing(opt, chunkFile, xCenters, yCenters);

    poreMask = validMask & (ncountGrid < opt.ThresholdN);
    matrixMask = validMask & ~poreMask;
    phaseGrid = nan(size(ncountGrid));
    phaseGrid(matrixMask) = 0;
    phaseGrid(poreMask) = 1;

    pore = analyzePhase('pore', poreMask, validMask, xCenters, yCenters, dx, dy, opt);
    matrix = analyzePhase('matrix', matrixMask, validMask, xCenters, yCenters, dx, dy, opt);
    profile = buildDirectionalProfiles(poreMask, validMask, xCenters, yCenters, dx, dy, opt, pore);

    cellArea = dx * dy;
    numValidCells = nnz(validMask);
    validArea = numValidCells * cellArea;
    poreCells = nnz(poreMask);
    matrixCells = nnz(matrixMask);
    poreArea = poreCells * cellArea;
    matrixArea = matrixCells * cellArea;
    interfaceLength = pore.totalInterfacePerimeter;

    globalStats = struct();
    globalStats.presentRows = presentCount;
    globalStats.numValidCells = numValidCells;
    globalStats.validArea = validArea;
    globalStats.poreCells = poreCells;
    globalStats.poreArea = poreArea;
    globalStats.matrixCells = matrixCells;
    globalStats.matrixArea = matrixArea;
    globalStats.porosity = safeDivide(poreArea, validArea);
    globalStats.matrixFraction = safeDivide(matrixArea, validArea);
    globalStats.interfaceLength = interfaceLength;
    globalStats.specificInterfaceBulk = safeDivide(interfaceLength, validArea);
    globalStats.specificInterfacePore = safeDivide(interfaceLength, poreArea);
    globalStats.specificInterfaceMatrix = safeDivide(interfaceLength, matrixArea);

    stats = buildPaperStats(globalStats, pore, matrix, poreMask, matrixMask, validMask, ...
        xCenters, yCenters, dx, dy, step.timestep, chunkFile, opt);

    plots = emptyPlotsStruct();
    if opt.MakePlots
        plots.phaseFig = plotPhaseGrid(xCenters, yCenters, phaseGrid, step.timestep, ...
            opt.ThresholdN, opt.PlotRangeX, opt.PlotRangeY);
        plots.poreLabelFig = plotLabelGrid(xCenters, yCenters, pore.labelGrid, validMask, ...
            sprintf('Pore Labels @ timestep %g', step.timestep), opt.PlotRangeX, opt.PlotRangeY);
        plots.matrixLabelFig = plotLabelGrid(xCenters, yCenters, matrix.labelGrid, validMask, ...
            sprintf('Matrix Labels @ timestep %g', step.timestep), opt.PlotRangeX, opt.PlotRangeY);
        plots.connectivityFig = plotConnectivityHighlights(xCenters, yCenters, validMask, ...
            pore, matrix, step.timestep, opt.PlotRangeX, opt.PlotRangeY);

        plots.poreCountFig = plotCountDistribution(pore.equivDiameterDistribution, 'Pore');
        plots.matrixCountFig = plotCountDistribution(matrix.equivDiameterDistribution, 'Matrix');

        if strcmp(opt.PositionAxis, 'x') || strcmp(opt.PositionAxis, 'both')
            plots.porePositionFigX = plotPositionDistribution(pore.positionDistribution.x, 'Pore', 'x');
            plots.matrixPositionFigX = plotPositionDistribution(matrix.positionDistribution.x, 'Matrix', 'x');
        end
        if strcmp(opt.PositionAxis, 'y') || strcmp(opt.PositionAxis, 'both')
            plots.porePositionFigY = plotPositionDistribution(pore.positionDistribution.y, 'Pore', 'y');
            plots.matrixPositionFigY = plotPositionDistribution(matrix.positionDistribution.y, 'Matrix', 'y');
        end
        if isfield(profile, 'x')
            plots.profileFigX = plotDirectionalProfile(profile.x, step.timestep);
        end
        if isfield(profile, 'y')
            plots.profileFigY = plotDirectionalProfile(profile.y, step.timestep);
        end
    end

    if opt.EnableEvolution
        stats.evolution = buildEvolutionSeries(chunkFile, opt);
        if opt.MakePlots
            plots.evolutionFig = plotEvolutionSeries(stats.evolution);
        end
    else
        stats.evolution = emptyEvolutionStats();
    end

    summary = struct();
    summary.porosity = globalStats.porosity;
    summary.specificInterfaceBulk = globalStats.specificInterfaceBulk;
    summary.poreNumComponents = pore.numComponents;
    summary.matrixNumComponents = matrix.numComponents;
    summary.largestPoreArea = pore.largestComponentArea;
    summary.largestMatrixArea = matrix.largestComponentArea;
    summary.poreConnectedX = pore.connectivity.x.isConnected;
    summary.poreConnectedY = pore.connectivity.y.isConnected;
    summary.matrixConnectedX = matrix.connectivity.x.isConnected;
    summary.matrixConnectedY = matrix.connectivity.y.isConnected;
    summary.profileAxis = opt.ProfileAxis;
    summary.profileNumBins = round(opt.ProfileNumBins);

    out = struct();
    out.filePath = chunkFile;
    out.selection = selectionInfo;
    out.stepIndex = step.stepIndex;
    out.timestep = step.timestep;
    out.thresholdN = opt.ThresholdN;
    out.connectivity = opt.Connectivity;
    out.boundary = opt.Boundary;
    out.coordScale = opt.CoordScale;
    out.ncountVar = opt.NcountVar;
    out.positionRange = struct('x', opt.PositionRangeX, 'y', opt.PositionRangeY);
    out.plotRange = struct('x', opt.PlotRangeX, 'y', opt.PlotRangeY);
    out.profileRange = struct('x', opt.ProfileRangeX, 'y', opt.ProfileRangeY);
    out.xCenters = xCenters;
    out.yCenters = yCenters;
    out.dx = dx;
    out.dy = dy;
    out.spacingSource = spacingSource;
    out.NcountGrid = ncountGrid;
    out.validMask = validMask;
    out.poreMask = poreMask;
    out.matrixMask = matrixMask;
    out.global = globalStats;
    out.summary = summary;
    out.profile = profile;
    out.pore = pore;
    out.matrix = matrix;
    out.stats = stats;
    out.plots = plots;
end

function validateInputs(opt)
    if isempty(opt.ThresholdN) || ~isscalar(opt.ThresholdN) || ~isfinite(opt.ThresholdN)
        error('analyze_chunk_network2d:MissingThresholdN', ...
            'ThresholdN must be provided as a finite scalar.');
    end
    if ~(isscalar(opt.Connectivity) && isfinite(opt.Connectivity) && round(opt.Connectivity) == 4)
        error('analyze_chunk_network2d:BadConnectivity', ...
            'Only Connectivity=4 is supported in this version.');
    end
    validBoundary = {'open', 'periodic-x', 'periodic-y', 'periodic-xy'};
    if ~any(strcmp(opt.Boundary, validBoundary))
        error('analyze_chunk_network2d:BadBoundary', ...
            'Boundary must be open/periodic-x/periodic-y/periodic-xy.');
    end
    validAxis = {'x', 'y', 'both'};
    if ~any(strcmp(opt.PositionAxis, validAxis))
        error('analyze_chunk_network2d:BadPositionAxis', ...
            'PositionAxis must be x/y/both.');
    end
    if ~any(strcmp(opt.ProfileAxis, validAxis))
        error('analyze_chunk_network2d:BadProfileAxis', ...
            'ProfileAxis must be x/y/both.');
    end
    if ~(isscalar(opt.PositionNumBins) && isfinite(opt.PositionNumBins) && opt.PositionNumBins >= 1)
        error('analyze_chunk_network2d:BadPositionNumBins', ...
            'PositionNumBins must be a positive scalar.');
    end
    if ~(isscalar(opt.ProfileNumBins) && isfinite(opt.ProfileNumBins) && opt.ProfileNumBins >= 1)
        error('analyze_chunk_network2d:BadProfileNumBins', ...
            'ProfileNumBins must be a positive scalar.');
    end
    if ~(isscalar(opt.HistNumBins) && isfinite(opt.HistNumBins) && opt.HistNumBins >= 1)
        error('analyze_chunk_network2d:BadHistNumBins', ...
            'HistNumBins must be a positive scalar.');
    end
    validateMeanPower(opt.MeanPowerM, 'MeanPowerM');
    validateMeanPower(opt.MeanPowerN, 'MeanPowerN');
    if ~(isscalar(opt.CoordScale) && isnumeric(opt.CoordScale) && isfinite(opt.CoordScale) && opt.CoordScale > 0)
        error('analyze_chunk_network2d:BadCoordScale', ...
            'CoordScale must be a positive finite scalar.');
    end
    validEvolutionAxis = {'', 'index', 'timestep', 'time'};
    if ~any(strcmp(opt.EvolutionSelectBy, validEvolutionAxis))
        error('analyze_chunk_network2d:BadEvolutionSelectBy', ...
            'EvolutionSelectBy must be empty or Index/TimeStep/Time.');
    end
    if ~(isscalar(opt.EvolutionStride) && isfinite(opt.EvolutionStride) && opt.EvolutionStride >= 1)
        error('analyze_chunk_network2d:BadEvolutionStride', ...
            'EvolutionStride must be a positive scalar.');
    end
    validatePositiveScalarOrEmpty(opt.Dx, 'Dx');
    validatePositiveScalarOrEmpty(opt.Dy, 'Dy');
    validateRangeOrEmpty(opt.PositionRangeX, 'PositionRangeX');
    validateRangeOrEmpty(opt.PositionRangeY, 'PositionRangeY');
    validateRangeOrEmpty(opt.PlotRangeX, 'PlotRangeX');
    validateRangeOrEmpty(opt.PlotRangeY, 'PlotRangeY');
    validateRangeOrEmpty(opt.ProfileRangeX, 'ProfileRangeX');
    validateRangeOrEmpty(opt.ProfileRangeY, 'ProfileRangeY');
    validateRangeOrEmpty(opt.EvolutionRange, 'EvolutionRange');
end

function validatePositiveScalarOrEmpty(v, name)
    if isempty(v)
        return;
    end
    if ~(isscalar(v) && isnumeric(v) && isfinite(v) && (v > 0))
        error('analyze_chunk_network2d:BadSpacing', ...
            '%s must be empty or a positive finite scalar.', name);
    end
end

function validateRangeOrEmpty(v, name)
    if isempty(v)
        return;
    end
    if ~(isnumeric(v) && numel(v) == 2 && all(isfinite(v(:))) && (v(2) >= v(1)))
        error('analyze_chunk_network2d:BadRange', ...
            '%s must be empty or a finite [min max] range with max >= min.', name);
    end
end

function validateMeanPower(v, name)
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v))
        error('analyze_chunk_network2d:BadMeanPower', ...
            '%s must be a finite real scalar.', name);
    end
end

function [x, y] = scaleCoordinates(x, y, coordScale)
    x = x .* coordScale;
    y = y .* coordScale;
end

function [x, y, ncount] = extractRequiredColumns(data, col, ncountVar)
    x = [];
    y = [];
    if isfield(col, 'Coord1')
        x = data(:, col.Coord1);
    elseif isfield(col, 'c_x')
        x = data(:, col.c_x);
    end

    if isfield(col, 'Coord2')
        y = data(:, col.Coord2);
    elseif isfield(col, 'c_y')
        y = data(:, col.c_y);
    end

    if isempty(x) || isempty(y)
        error('analyze_chunk_network2d:Not2DChunk', ...
            '2D chunk data requires x/y coordinates (Coord1/Coord2 or c_x/c_y).');
    end
    if ~isfield(col, ncountVar)
        error('analyze_chunk_network2d:MissingNcount', ...
            'Missing Ncount variable "%s".', ncountVar);
    end

    ncount = data(:, col.(ncountVar));
    if any(~isfinite(x)) || any(~isfinite(y)) || any(~isfinite(ncount))
        error('analyze_chunk_network2d:NonFiniteInput', ...
            'x, y, and Ncount must be finite for every present chunk cell.');
    end
end

function [ncountGrid, validMask, xCenters, yCenters, presentCount] = buildGrid(x, y, ncount)
    xCenters = unique(x(:)).';
    yCenters = unique(y(:)).';
    nx = numel(xCenters);
    ny = numel(yCenters);

    [~, ix] = ismember(x(:), xCenters(:));
    [~, iy] = ismember(y(:), yCenters(:));
    lin = sub2ind([ny, nx], iy, ix);

    if numel(unique(lin)) ~= numel(lin)
        error('analyze_chunk_network2d:DuplicateGridCell', ...
            'Duplicate coordinate pairs detected in the 2D chunk data.');
    end

    ncountGrid = nan(ny, nx);
    validMask = false(ny, nx);
    ncountGrid(lin) = ncount(:);
    validMask(lin) = true;
    presentCount = numel(lin);
end

function [dx, dy, source] = resolveGridSpacing(opt, chunkFile, xCenters, yCenters)
    meta = parseFilenameMetadata(chunkFile);
    [dx, srcDx] = resolveOneSpacing(opt.Dx, 'dx', xCenters, meta, opt.CoordScale);
    [dy, srcDy] = resolveOneSpacing(opt.Dy, 'dy', yCenters, meta, opt.CoordScale);
    source = struct('dx', srcDx, 'dy', srcDy);
end

function [spacing, src] = resolveOneSpacing(userValue, key, coords, meta, coordScale)
    if ~isempty(userValue)
        spacing = userValue .* coordScale;
        src = 'option';
        return;
    end
    if isfield(meta, key)
        spacing = meta.(key) * coordScale;
        if isfinite(spacing) && spacing > 0
            src = 'filename';
            return;
        end
    end
    spacing = inferSpacingFromCoords(coords, key);
    src = 'coordinates';
end

function spacing = inferSpacingFromCoords(coords, key)
    coords = coords(:);
    if numel(coords) < 2
        error('analyze_chunk_network2d:NeedSpacing', ...
            'Cannot infer %s from coordinates with fewer than 2 unique centers.', key);
    end

    diffs = diff(sort(coords));
    if any(~isfinite(diffs)) || any(diffs <= 0)
        error('analyze_chunk_network2d:BadCoordinateSpacing', ...
            'Cannot infer %s from non-monotone coordinates.', key);
    end

    spacing = median(diffs);
    tol = max(1e-9, 1e-2 * max(abs(spacing), 1));
    if any(abs(diffs - spacing) > tol)
        error('analyze_chunk_network2d:IrregularGrid', ...
            'Cannot infer %s from irregular coordinate spacing.', key);
    end
end

function meta = parseFilenameMetadata(chunkFile)
    [~, nameOnly, ~] = fileparts(chunkFile);
    parts = strsplit(nameOnly, '_');
    meta = struct();
    if numel(parts) < 2
        return;
    end

    for i = 2:(numel(parts) - 1)
        key = matlab.lang.makeValidName(parts{i});
        val = str2double(parts{i+1});
        if ~isnan(val)
            meta.(key) = val;
        end
    end
end

function phase = analyzePhase(name, phaseMask, validMask, xCenters, yCenters, dx, dy, opt)
    [labelGrid, wrapXComp, wrapYComp] = labelConnectedComponents(phaseMask, validMask, opt.Boundary);
    numComponents = max(labelGrid(:));
    comp = computeComponentStats(labelGrid, numComponents, phaseMask, validMask, ...
        xCenters, yCenters, dx, dy, opt.Boundary, wrapXComp, wrapYComp);

    area = nnz(phaseMask) * dx * dy;
    validArea = nnz(validMask) * dx * dy;
    largestArea = 0;
    if ~isempty(comp.area)
        largestArea = max(comp.area);
    end

    hasMissingCells = any(~validMask(:));
    topology = computePhaseTopologyDigital(phaseMask, opt.Boundary, numComponents, hasMissingCells);
    holeCount = topology.beta1;
    eulerCharacteristic = topology.chi;
    topologyNote = topology.note;

    areaStats = computeSummaryStats(comp.area);
    diamStats = computeSummaryStats(comp.equivDiameter, opt.MeanPowerM, opt.MeanPowerN);
    equivDiameterDistribution = buildHistogramDistribution(comp.equivDiameter, round(opt.HistNumBins));
    positionDistribution = struct();
    positionDistribution.x = buildPositionDistribution(comp.centroidX, comp.equivDiameter, ...
        round(opt.PositionNumBins), opt.PositionRangeX, dx, opt.MeanPowerM, opt.MeanPowerN);
    positionDistribution.y = buildPositionDistribution(comp.centroidY, comp.equivDiameter, ...
        round(opt.PositionNumBins), opt.PositionRangeY, dy, opt.MeanPowerM, opt.MeanPowerN);
    connectivity = buildDirectionalConnectivity(comp, opt.Boundary, area);

    phase = struct();
    phase.name = name;
    phase.mask = phaseMask;
    phase.labelGrid = labelGrid;
    phase.numCells = nnz(phaseMask);
    phase.area = area;
    phase.areaFraction = safeDivide(area, validArea);
    phase.numComponents = numComponents;
    phase.largestComponentArea = largestArea;
    phase.largestComponentFraction = safeDivideOrZero(largestArea, area);
    phase.percolatesX = computePercolates(comp, 'x', opt.Boundary);
    phase.percolatesY = computePercolates(comp, 'y', opt.Boundary);
    phase.wrapsX = computeWraps(comp, 'x', opt.Boundary);
    phase.wrapsY = computeWraps(comp, 'y', opt.Boundary);
    phase.holeCount = holeCount;
    phase.eulerCharacteristic = eulerCharacteristic;
    phase.topologyBeta0 = topology.beta0;
    phase.topologyBeta1 = topology.beta1;
    phase.topologyChi = topology.chi;
    phase.topologyBeta2 = topology.beta2;
    phase.topologyNote = topologyNote;
    phase.totalPerimeterOpen = sum(comp.perimeterOpen);
    phase.totalInterfacePerimeter = sum(comp.interfacePerimeter);
    phase.components = comp;
    phase.areaStats = areaStats;
    phase.equivDiameterStats = diamStats;
    phase.equivDiameterDistribution = equivDiameterDistribution;
    phase.positionDistribution = positionDistribution;
    phase.connectivity = connectivity;
end

function [labelGrid, wrapXComp, wrapYComp] = labelConnectedComponents(phaseMask, validMask, boundary)
    [ny, nx] = size(phaseMask);
    labelGrid = zeros(ny, nx);
    numPhase = nnz(phaseMask);
    if numPhase == 0
        wrapXComp = false(0, 1);
        wrapYComp = false(0, 1);
        return;
    end

    queue = zeros(numPhase, 1);
    wrapXComp = false(numPhase, 1);
    wrapYComp = false(numPhase, 1);
    offsetX = nan(ny, nx);
    offsetY = nan(ny, nx);
    compId = 0;
    phaseLin = find(phaseMask);

    for idx = 1:numel(phaseLin)
        seed = phaseLin(idx);
        if labelGrid(seed) ~= 0
            continue;
        end

        compId = compId + 1;
        head = 1;
        tail = 1;
        queue(1) = seed;
        labelGrid(seed) = compId;
        offsetX(seed) = 0;
        offsetY(seed) = 0;
        wrapX = false;
        wrapY = false;

        while head <= tail
            cur = queue(head);
            head = head + 1;
            [r, c] = ind2sub([ny, nx], cur);
            curOffsetX = offsetX(cur);
            curOffsetY = offsetY(cur);

            for d = 1:4
                [nr, nc, hasNeighbor, shiftX, shiftY] = neighborAtWithShift(r, c, d, ny, nx, boundary);
                if ~hasNeighbor
                    continue;
                end
                if ~validMask(nr, nc) || ~phaseMask(nr, nc)
                    continue;
                end

                nextOffsetX = curOffsetX + shiftX;
                nextOffsetY = curOffsetY + shiftY;

                if labelGrid(nr, nc) == 0
                    tail = tail + 1;
                    nextLin = sub2ind([ny, nx], nr, nc);
                    queue(tail) = nextLin;
                    labelGrid(nextLin) = compId;
                    offsetX(nextLin) = nextOffsetX;
                    offsetY(nextLin) = nextOffsetY;
                else
                    deltaX = nextOffsetX - offsetX(nr, nc);
                    deltaY = nextOffsetY - offsetY(nr, nc);
                    wrapX = wrapX || (deltaX ~= 0);
                    wrapY = wrapY || (deltaY ~= 0);
                end
            end
        end

        wrapXComp(compId) = wrapX;
        wrapYComp(compId) = wrapY;
    end

    wrapXComp = wrapXComp(1:compId);
    wrapYComp = wrapYComp(1:compId);
end

function comp = computeComponentStats(labelGrid, numComponents, phaseMask, validMask, ...
        xCenters, yCenters, dx, dy, boundary, wrapXComp, wrapYComp)
    comp = emptyComponentStats();
    if numComponents == 0
        return;
    end

    [ny, nx] = size(labelGrid);
    comp.label = (1:numComponents).';
    comp.cellCount = zeros(numComponents, 1);
    comp.area = zeros(numComponents, 1);
    comp.equivDiameter = zeros(numComponents, 1);
    comp.centroidX = zeros(numComponents, 1);
    comp.centroidY = zeros(numComponents, 1);
    comp.bboxWidth = zeros(numComponents, 1);
    comp.bboxHeight = zeros(numComponents, 1);
    comp.perimeterOpen = zeros(numComponents, 1);
    comp.interfacePerimeter = zeros(numComponents, 1);
    comp.shapeFactor = nan(numComponents, 1);
    comp.touchesLeft = false(numComponents, 1);
    comp.touchesRight = false(numComponents, 1);
    comp.touchesBottom = false(numComponents, 1);
    comp.touchesTop = false(numComponents, 1);
    comp.wrapX = false(numComponents, 1);
    comp.wrapY = false(numComponents, 1);
    comp.percolatesX = false(numComponents, 1);
    comp.percolatesY = false(numComponents, 1);

    for k = 1:numComponents
        mask = (labelGrid == k);
        [rows, cols] = find(mask);
        comp.cellCount(k) = numel(rows);
        comp.area(k) = comp.cellCount(k) * dx * dy;
        comp.equivDiameter(k) = 2 * sqrt(comp.area(k) / pi);
        comp.centroidX(k) = mean(xCenters(cols));
        comp.centroidY(k) = mean(yCenters(rows));
        comp.bboxWidth(k) = (max(cols) - min(cols) + 1) * dx;
        comp.bboxHeight(k) = (max(rows) - min(rows) + 1) * dy;

        comp.touchesLeft(k) = any(cols == 1);
        comp.touchesRight(k) = any(cols == nx);
        comp.touchesBottom(k) = any(rows == 1);
        comp.touchesTop(k) = any(rows == ny);
        comp.wrapX(k) = wrapXComp(k);
        comp.wrapY(k) = wrapYComp(k);
        comp.percolatesX(k) = comp.touchesLeft(k) && comp.touchesRight(k);
        comp.percolatesY(k) = comp.touchesBottom(k) && comp.touchesTop(k);

        [comp.perimeterOpen(k), comp.interfacePerimeter(k)] = ...
            computePerimeterForMask(mask, phaseMask, validMask, dx, dy, boundary);
        if comp.perimeterOpen(k) > 0
            comp.shapeFactor(k) = 4 * pi * comp.area(k) / (comp.perimeterOpen(k) ^ 2);
        end
    end
end

function [perimeterOpen, interfacePerimeter] = computePerimeterForMask(componentMask, phaseMask, validMask, dx, dy, boundary)
    [ny, nx] = size(componentMask);
    perimeterOpen = 0;
    interfacePerimeter = 0;
    [rows, cols] = find(componentMask);

    for i = 1:numel(rows)
        r = rows(i);
        c = cols(i);
        for d = 1:4
            [nr, nc, hasNeighbor, ~, ~] = neighborAt(r, c, d, ny, nx, boundary);
            sideLength = sideLengthByDir(d, dx, dy);
            if ~hasNeighbor
                perimeterOpen = perimeterOpen + sideLength;
                continue;
            end
            if ~validMask(nr, nc)
                perimeterOpen = perimeterOpen + sideLength;
                continue;
            end
            if phaseMask(nr, nc)
                continue;
            end

            perimeterOpen = perimeterOpen + sideLength;
            interfacePerimeter = interfacePerimeter + sideLength;
        end
    end
end

function value = computePercolates(comp, axisName, boundary)
    if strcmp(axisName, 'x')
        if hasPeriodicX(boundary)
            value = NaN;
        else
            value = any(comp.percolatesX);
        end
    else
        if hasPeriodicY(boundary)
            value = NaN;
        else
            value = any(comp.percolatesY);
        end
    end
end

function value = computeWraps(comp, axisName, boundary)
    if strcmp(axisName, 'x')
        if hasPeriodicX(boundary)
            value = any(comp.wrapX);
        else
            value = false;
        end
    else
        if hasPeriodicY(boundary)
            value = any(comp.wrapY);
        else
            value = false;
        end
    end
end

function profiles = buildDirectionalProfiles(poreMask, validMask, xCenters, yCenters, dx, dy, opt, porePhase)
    profiles = struct();
    if strcmp(opt.ProfileAxis, 'x') || strcmp(opt.ProfileAxis, 'both')
        profiles.x = buildDirectionalProfile(poreMask, validMask, porePhase, ...
            xCenters, yCenters, dx, dy, 'x', round(opt.ProfileNumBins), opt);
    end
    if strcmp(opt.ProfileAxis, 'y') || strcmp(opt.ProfileAxis, 'both')
        profiles.y = buildDirectionalProfile(poreMask, validMask, porePhase, ...
            xCenters, yCenters, dx, dy, 'y', round(opt.ProfileNumBins), opt);
    end
end

function profile = buildDirectionalProfile(poreMask, validMask, porePhase, ...
        xCenters, yCenters, dx, dy, axisName, nBins, opt)
    cellArea = dx * dy;
    [axisCenters, coordMin, coordMax, perpAxis, axisStep] = resolveProfileAxis(axisName, xCenters, yCenters, dx, dy);
    coordRange = getProfileRange(opt, axisName);
    edges = buildAxisEdges(coordRange, coordMin, coordMax, axisStep, nBins);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    binIdCenter = discretize(axisCenters(:), edges);
    interfaceLength = computeInterfaceLengthProfile(poreMask, validMask, xCenters, yCenters, dx, dy, axisName, edges);

    validArea = zeros(1, nBins);
    poreArea = zeros(1, nBins);
    matrixArea = zeros(1, nBins);
    porosity = nan(1, nBins);
    specificInterfaceBulk = nan(1, nBins);
    specificInterfacePore = nan(1, nBins);
    specificInterfaceMatrix = nan(1, nBins);

    for i = 1:nBins
        idxAxis = find(binIdCenter == i);
        if isempty(idxAxis)
            continue;
        end

        [poreSub, validSub] = extractSliceSubgrid(poreMask, validMask, axisName, idxAxis);
        validArea(i) = nnz(validSub) * cellArea;
        poreArea(i) = nnz(poreSub) * cellArea;
        matrixArea(i) = nnz(validSub & ~poreSub) * cellArea;
        porosity(i) = safeDivide(poreArea(i), validArea(i));
        specificInterfaceBulk(i) = safeDivide(interfaceLength(i), validArea(i));
        specificInterfacePore(i) = safeDivide(interfaceLength(i), poreArea(i));
        specificInterfaceMatrix(i) = safeDivide(interfaceLength(i), matrixArea(i));

    end

    connProfile = buildComponentConnectivityProfile(porePhase.components, porePhase.connectivity, ...
        axisName, perpAxis, edges);

    profile = struct();
    profile.axis = axisName;
    profile.perpendicularAxis = perpAxis;
    profile.range = [edges(1), edges(end)];
    profile.edges = edges;
    profile.centers = centers;
    profile.validArea = validArea;
    profile.poreArea = poreArea;
    profile.matrixArea = matrixArea;
    profile.porosity = porosity;
    profile.interfaceLength = interfaceLength;
    profile.specificInterfaceBulk = specificInterfaceBulk;
    profile.specificInterfacePore = specificInterfacePore;
    profile.specificInterfaceMatrix = specificInterfaceMatrix;
    profile.connectivityFlag = connProfile.flag;
    profile.connectivityFraction = connProfile.fraction;
    profile.connectivityArea = connProfile.area;
    profile.connectivityComponentCount = connProfile.count;
    profile.componentPoreArea = connProfile.componentPoreArea;
end

function [axisCenters, coordMin, coordMax, perpAxis, axisStep] = resolveProfileAxis(axisName, xCenters, yCenters, dx, dy)
    if strcmp(axisName, 'x')
        axisCenters = xCenters;
        coordMin = xCenters(1) - 0.5 * dx;
        coordMax = xCenters(end) + 0.5 * dx;
        perpAxis = 'y';
        axisStep = dx;
    else
        axisCenters = yCenters;
        coordMin = yCenters(1) - 0.5 * dy;
        coordMax = yCenters(end) + 0.5 * dy;
        perpAxis = 'x';
        axisStep = dy;
    end
end

function interfaceLength = computeInterfaceLengthProfile(poreMask, validMask, xCenters, yCenters, dx, dy, axisName, edges)
    [ny, nx] = size(poreMask);
    nBins = numel(edges) - 1;
    interfaceLength = zeros(1, nBins);

    for r = 1:ny-1
        for c = 1:nx
            if ~(validMask(r, c) && validMask(r+1, c))
                continue;
            end
            if poreMask(r, c) == poreMask(r+1, c)
                continue;
            end
            if strcmp(axisName, 'x')
                coord = xCenters(c);
            else
                coord = 0.5 * (yCenters(r) + yCenters(r+1));
            end
            binId = locateBin(coord, edges);
            if binId > 0
                interfaceLength(binId) = interfaceLength(binId) + dx;
            end
        end
    end

    for r = 1:ny
        for c = 1:nx-1
            if ~(validMask(r, c) && validMask(r, c+1))
                continue;
            end
            if poreMask(r, c) == poreMask(r, c+1)
                continue;
            end
            if strcmp(axisName, 'x')
                coord = 0.5 * (xCenters(c) + xCenters(c+1));
            else
                coord = yCenters(r);
            end
            binId = locateBin(coord, edges);
            if binId > 0
                interfaceLength(binId) = interfaceLength(binId) + dy;
            end
        end
    end
end

function [poreSub, validSub] = extractSliceSubgrid(poreMask, validMask, axisName, idxAxis)
    if strcmp(axisName, 'x')
        poreSub = poreMask(:, idxAxis);
        validSub = validMask(:, idxAxis);
    else
        poreSub = poreMask(idxAxis, :);
        validSub = validMask(idxAxis, :);
    end
end

function range = getProfileRange(opt, axisName)
    if strcmp(axisName, 'x')
        range = opt.ProfileRangeX;
    else
        range = opt.ProfileRangeY;
    end
end

function edges = buildAxisEdges(coordRange, coordMin, coordMax, axisStep, nBins)
    if isempty(coordRange)
        minEdge = coordMin;
        maxEdge = coordMax;
    else
        minEdge = coordRange(1);
        maxEdge = coordRange(2);
    end
    if minEdge == maxEdge
        pad = 0.5 * axisStep;
        minEdge = minEdge - pad;
        maxEdge = maxEdge + pad;
    end
    edges = linspace(minEdge, maxEdge, nBins + 1);
end

function profile = buildComponentConnectivityProfile(comp, connectivity, axisName, perpAxis, edges)
    nBins = numel(edges) - 1;
    profile = struct('flag', false(1, nBins), ...
        'fraction', zeros(1, nBins), ...
        'area', zeros(1, nBins), ...
        'count', zeros(1, nBins), ...
        'componentPoreArea', zeros(1, nBins));
    if isempty(comp.label)
        return;
    end

    if strcmp(axisName, 'x')
        centroid = comp.centroidX;
    else
        centroid = comp.centroidY;
    end
    binId = nan(size(centroid));
    for k = 1:numel(centroid)
        binId(k) = locateBin(centroid(k), edges);
    end
    connectedIds = getConnectivityIds(connectivity, perpAxis);
    connectedMask = ismember(comp.label, connectedIds);

    for i = 1:nBins
        idxBin = (binId == i);
        if ~any(idxBin)
            continue;
        end
        profile.componentPoreArea(i) = sum(comp.area(idxBin));
        idxConnected = idxBin & connectedMask;
        profile.flag(i) = any(idxConnected);
        profile.count(i) = sum(idxConnected);
        profile.area(i) = sum(comp.area(idxConnected));
        profile.fraction(i) = safeDivide(profile.area(i), profile.componentPoreArea(i));
        if isnan(profile.fraction(i))
            profile.fraction(i) = 0;
        end
    end
end

function ids = getConnectivityIds(connectivity, axisName)
    if strcmp(axisName, 'x')
        ids = connectivity.x.componentIds;
    else
        ids = connectivity.y.componentIds;
    end
end

function binId = locateBin(coord, edges)
    binId = discretize(coord, edges);
    if isnan(binId)
        tol = max(1e-12, 1e-12 * max(abs(edges(end)), 1));
        if abs(coord - edges(end)) <= tol
            binId = numel(edges) - 1;
        else
            binId = 0;
        end
    end
end

function connectivity = buildDirectionalConnectivity(comp, boundary, phaseArea)
    connectivity = struct();
    connectivity.x = buildConnectivityAxis(comp, boundary, phaseArea, 'x');
    connectivity.y = buildConnectivityAxis(comp, boundary, phaseArea, 'y');
end

function info = buildConnectivityAxis(comp, boundary, phaseArea, axisName)
    if strcmp(axisName, 'x')
        if hasPeriodicX(boundary)
            componentMask = comp.wrapX;
            criterion = 'wraps-periodic-boundary';
        else
            componentMask = comp.percolatesX;
            criterion = 'touches-both-open-boundaries';
        end
    else
        if hasPeriodicY(boundary)
            componentMask = comp.wrapY;
            criterion = 'wraps-periodic-boundary';
        else
            componentMask = comp.percolatesY;
            criterion = 'touches-both-open-boundaries';
        end
    end

    ids = comp.label(componentMask);
    areas = comp.area(componentMask);

    info = struct();
    info.axis = axisName;
    info.criterion = criterion;
    info.isConnected = any(componentMask);
    info.componentIds = ids(:).';
    info.componentCount = numel(ids);
    info.totalConnectedArea = sum(areas);
    info.connectedAreaFraction = safeDivideOrZero(sum(areas), phaseArea);
    if isempty(areas)
        info.largestConnectedArea = 0;
    else
        info.largestConnectedArea = max(areas);
    end
end

function stats = computeSummaryStats(x, meanPowerM, meanPowerN)
    x = x(:);
    x = x(isfinite(x));

    if nargin < 2
        meanPowerM = 1;
        meanPowerN = 0;
    end

    stats = struct();
    stats.n = numel(x);
    if isempty(x)
        stats.min = NaN;
        stats.max = NaN;
        stats.mean = NaN;
        stats.std = NaN;
        stats.median = NaN;
        stats.quantileProb = [0.05, 0.25, 0.50, 0.75, 0.95];
        stats.quantileValue = nan(1, 5);
        return;
    end

    stats.min = min(x);
    stats.max = max(x);
    stats.mean = momentRatioMean(x, meanPowerM, meanPowerN);
    stats.std = std(x);
    stats.median = median(x);
    stats.quantileProb = [0.05, 0.25, 0.50, 0.75, 0.95];
    stats.quantileValue = computeQuantiles(x, stats.quantileProb);
end

function q = computeQuantiles(x, prob)
    x = sort(x(:), 'ascend');
    n = numel(x);
    q = nan(size(prob));
    if n == 0
        return;
    end
    if n == 1
        q(:) = x(1);
        return;
    end

    for i = 1:numel(prob)
        p = min(max(prob(i), 0), 1);
        pos = 1 + (n - 1) * p;
        lo = floor(pos);
        hi = ceil(pos);
        if lo == hi
            q(i) = x(lo);
        else
            t = pos - lo;
            q(i) = x(lo) * (1 - t) + x(hi) * t;
        end
    end
end

function dist = buildHistogramDistribution(x, nBins)
    x = x(:);
    x = x(isfinite(x) & (x > 0));
    dist = struct('edges', [], 'centers', [], 'count', [], 'probability', [], ...
        'cdfX', [], 'cdfCount', [], 'cdfProbability', []);

    if isempty(x)
        return;
    end

    nBins = max(1, round(nBins));
    xmin = min(x);
    xmax = max(x);
    if xmin == xmax
        pad = max(1e-12, abs(xmin) * 1e-6);
        edges = linspace(xmin - pad, xmax + pad, nBins + 1);
    else
        edges = linspace(xmin, xmax, nBins + 1);
    end

    count = histcounts(x, edges);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    probability = count ./ sum(count);

    xs = sort(x, 'ascend');
    cdfCount = (1:numel(xs)).';
    cdfProbability = cdfCount ./ numel(xs);

    dist.edges = edges;
    dist.centers = centers;
    dist.count = count;
    dist.probability = probability;
    dist.cdfX = xs;
    dist.cdfCount = cdfCount;
    dist.cdfProbability = cdfProbability;
end

function out = buildPositionDistribution(position, sizeValue, nBins, coordRange, axisStep, meanPowerM, meanPowerN)
    position = position(:);
    sizeValue = sizeValue(:);
    if nargin < 6
        meanPowerM = 1;
        meanPowerN = 0;
    end
    valid = isfinite(position) & isfinite(sizeValue) & (sizeValue > 0);
    position = position(valid);
    sizeValue = sizeValue(valid);

    out = struct('edges', [], 'centers', [], 'meanEquivDiameter', [], 'count', []);
    if isempty(position)
        return;
    end

    if ~isempty(coordRange)
        inRange = (position >= coordRange(1)) & (position <= coordRange(2));
        position = position(inRange);
        sizeValue = sizeValue(inRange);
        if isempty(position)
            out.edges = buildAxisEdges(coordRange, coordRange(1), coordRange(2), axisStep, nBins);
            out.centers = 0.5 * (out.edges(1:end-1) + out.edges(2:end));
            out.meanEquivDiameter = nan(size(out.centers));
            out.count = zeros(size(out.centers));
            return;
        end
    end

    nBins = max(1, round(nBins));
    if isempty(coordRange)
        xmin = min(position);
        xmax = max(position);
        if xmin == xmax
            edges = [xmin - 0.5 * axisStep, xmax + 0.5 * axisStep];
        else
            edges = linspace(xmin, xmax, nBins + 1);
        end
    else
        edges = buildAxisEdges(coordRange, coordRange(1), coordRange(2), axisStep, nBins);
    end

    binId = nan(size(position));
    for i = 1:numel(position)
        binId(i) = locateBin(position(i), edges);
    end
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    meanSize = nan(size(centers));
    count = zeros(size(centers));
    for i = 1:numel(centers)
        idx = (binId == i);
        count(i) = sum(idx);
        if count(i) > 0
            meanSize(i) = momentRatioMean(sizeValue(idx), meanPowerM, meanPowerN);
        end
    end

    out.edges = edges;
    out.centers = centers;
    out.meanEquivDiameter = meanSize;
    out.count = count;
end

function meanVal = momentRatioMean(x, meanPowerM, meanPowerN)
    if isempty(x)
        meanVal = NaN;
        return;
    end

    num = sum(x .^ meanPowerM);
    den = sum(x .^ meanPowerN);
    if ~isfinite(num) || ~isfinite(den) || den == 0
        meanVal = NaN;
    else
        meanVal = num / den;
    end
end

function stats = buildPaperStats(globalStats, pore, matrix, poreMask, matrixMask, validMask, ...
        xCenters, yCenters, dx, dy, timestep, chunkFile, opt)
    toolbox = detectImageToolboxAvailability();

    stats = struct();
    stats.meta = buildPaperMeta(opt, dx, dy, timestep, chunkFile, xCenters, yCenters, toolbox, validMask);
    stats.topology = struct( ...
        'pore', buildTopologyPhaseStats(pore, opt.Boundary), ...
        'matrix', buildTopologyPhaseStats(matrix, opt.Boundary));
    stats.connectivity = struct( ...
        'pore', buildConnectivityPhaseStats(pore, opt.Boundary), ...
        'matrix', buildConnectivityPhaseStats(matrix, opt.Boundary));
    stats.geometry = buildGeometryStats(globalStats);
    stats.size = struct( ...
        'pore', buildPhaseSizeStats(pore, opt), ...
        'matrix', buildPhaseSizeStats(matrix, opt));
    stats.network = struct( ...
        'pore', buildPhaseSkeletonStats(poreMask, pore.area, dx, dy, opt.Boundary, opt.EnableSkeletonGraph, toolbox), ...
        'matrix', buildPhaseSkeletonStats(matrixMask, matrix.area, dx, dy, opt.Boundary, opt.EnableSkeletonGraph, toolbox));
    stats.thickness = struct( ...
        'matrix', buildThicknessStats(matrixMask, dx, dy, opt.Boundary, toolbox));
    stats.fragmentation = struct( ...
        'matrix', buildFragmentationStats(matrix, stats.size.matrix, stats.topology.matrix));
end

function meta = buildPaperMeta(opt, dx, dy, timestep, chunkFile, xCenters, yCenters, toolbox, validMask)
    meta = struct();
    meta.filePath = chunkFile;
    meta.timestep = timestep;
    meta.boundary = opt.Boundary;
    meta.thresholdN = opt.ThresholdN;
    meta.foregroundConnectivity = 4;
    meta.backgroundHoleConnectivity = 8;
    meta.meanDefinition = struct('m', opt.MeanPowerM, 'n', opt.MeanPowerN);
    meta.coordScale = opt.CoordScale;
    meta.dx = dx;
    meta.dy = dy;
    meta.gridSize = [numel(yCenters), numel(xCenters)];
    meta.hasMissingCells = any(~validMask(:));
    meta.missingCellCount = nnz(~validMask);
    meta.imageToolbox = toolbox;
    meta.topologyModel = ['Digital 2D cell complex consistent with 4-neighbor foreground: ', ...
        '0-cells are occupied pixels, 1-cells are occupied horizontal/vertical adjacencies, ', ...
        'and 2-cells are occupied 2x2 plaquettes.'];
    meta.periodicTopologyNote = ['Under periodic boundaries, beta1 counts all non-trivial 1D loops ', ...
        'of the wrapped cell complex, including ordinary enclosed holes and wrapping loops.'];
    if meta.hasMissingCells
        meta.missingCellNote = ['Missing chunk cells are treated as absent observed cells. ', ...
            'Loops around data gaps are therefore loops of the observed topology rather than inferred material.'];
    else
        meta.missingCellNote = '';
    end
end

function topo = buildTopologyPhaseStats(phase, boundary)
    topo = struct();
    topo.beta0 = phase.topologyBeta0;
    topo.beta1 = phase.topologyBeta1;
    topo.chi = phase.topologyChi;
    topo.beta2 = phase.topologyBeta2;
    topo.note = phase.topologyNote;
    if strcmp(boundary, 'open') && isempty(topo.note)
        topo.note = 'Under open boundaries, beta1 equals the ordinary hole count.';
    end
end

function topo = computePhaseTopologyDigital(phaseMask, boundary, beta0, hasMissingCells)
    phaseMask = logical(phaseMask);
    topo = struct();
    topo.beta0 = beta0;
    topo.beta1 = 0;
    topo.beta2 = 0;
    topo.chi = 0;
    topo.note = '';

    if ~any(phaseMask(:))
        topo.note = buildTopologyConventionNote(boundary, hasMissingCells);
        return;
    end

    n0 = nnz(phaseMask);
    n1 = countPhaseAdjacencyEdges(phaseMask, boundary);
    n2 = countPhasePlaquettes(phaseMask, boundary);
    chi = n0 - n1 + n2;
    beta2 = computePhaseBeta2(phaseMask, boundary);
    beta1 = beta0 + beta2 - chi;
    if beta1 < 0
        beta1 = 0;
    end

    topo.beta1 = beta1;
    topo.beta2 = beta2;
    topo.chi = chi;
    topo.note = buildTopologyConventionNote(boundary, hasMissingCells);
end

function nEdges = countPhaseAdjacencyEdges(phaseMask, boundary)
    [ny, nx] = size(phaseMask);
    nEdges = 0;

    for r = 1:ny
        for c = 1:nx
            if ~phaseMask(r, c)
                continue;
            end

            [nr, nc, hasNeighbor, ~, ~] = neighborAt(r, c, 2, ny, nx, boundary);
            if hasNeighbor && phaseMask(nr, nc)
                nEdges = nEdges + 1;
            end

            [nr, nc, hasNeighbor, ~, ~] = neighborAt(r, c, 4, ny, nx, boundary);
            if hasNeighbor && phaseMask(nr, nc)
                nEdges = nEdges + 1;
            end
        end
    end
end

function nFaces = countPhasePlaquettes(phaseMask, boundary)
    [ny, nx] = size(phaseMask);
    nFaces = 0;

    for r = 1:ny
        for c = 1:nx
            if ~phaseMask(r, c)
                continue;
            end

            [rr, ccRight, hasRight, ~, ~] = neighborAt(r, c, 2, ny, nx, boundary);
            [rrDown, cc, hasDown, ~, ~] = neighborAt(r, c, 4, ny, nx, boundary);
            if ~hasRight || ~hasDown
                continue;
            end
            if phaseMask(rr, ccRight) && phaseMask(rrDown, cc) && phaseMask(rrDown, ccRight)
                nFaces = nFaces + 1;
            end
        end
    end
end

function beta2 = computePhaseBeta2(phaseMask, boundary)
    beta2 = 0;
    if strcmp(boundary, 'periodic-xy') && all(phaseMask(:))
        beta2 = 1;
    end
end

function note = buildTopologyConventionNote(boundary, hasMissingCells)
    if strcmp(boundary, 'open')
        note = ['Topology uses the 4-neighbor digital cell complex. ', ...
            'Under open boundaries beta1 equals the ordinary hole count.'];
    else
        note = ['Topology uses the wrapped 4-neighbor digital cell complex. ', ...
            'Under periodic boundaries beta1 counts all non-trivial 1D loops, ', ...
            'including ordinary enclosed holes and periodic wrapping loops.'];
    end
    if hasMissingCells
        note = [note, ' Missing chunk cells are treated as absent observed cells and may create loops around data gaps.'];
    end
end

function conn = buildConnectivityPhaseStats(phase, boundary)
    conn = struct();
    conn.componentCount = phase.numComponents;
    conn.largestArea = phase.largestComponentArea;
    conn.largestFraction = phase.largestComponentFraction;
    conn.wrapsX = phase.wrapsX;
    conn.wrapsY = phase.wrapsY;
    if hasPeriodicX(boundary)
        conn.percolatesX = logical(phase.connectivity.x.isConnected);
    else
        conn.percolatesX = logical(phase.percolatesX);
    end
    if hasPeriodicY(boundary)
        conn.percolatesY = logical(phase.connectivity.y.isConnected);
    else
        conn.percolatesY = logical(phase.percolatesY);
    end
    conn.connectedAreaFractionX = phase.connectivity.x.connectedAreaFraction;
    conn.connectedAreaFractionY = phase.connectivity.y.connectedAreaFraction;
end

function geom = buildGeometryStats(globalStats)
    geom = struct();
    geom.phi = globalStats.porosity;
    geom.matrixFraction = globalStats.matrixFraction;
    geom.interfaceLength = globalStats.interfaceLength;
    geom.specificInterface = globalStats.specificInterfaceBulk;
    geom.validArea = globalStats.validArea;
    geom.poreArea = globalStats.poreArea;
    geom.matrixArea = globalStats.matrixArea;
end

function sizeStats = buildPhaseSizeStats(phase, opt)
    area = phase.components.area(:);
    diameter = phase.components.equivDiameter(:);
    radius = 0.5 .* diameter;

    sizeStats = struct();
    sizeStats.count = numel(area);
    sizeStats.area = area;
    sizeStats.radius = radius;
    sizeStats.diameter = diameter;
    sizeStats.maxArea = zeroIfEmpty(maxOrNaN(area));
    sizeStats.maxRadius = zeroIfEmpty(maxOrNaN(radius));
    sizeStats.maxDiameter = zeroIfEmpty(maxOrNaN(diameter));
    sizeStats.meanArea = momentRatioMean(area, opt.MeanPowerM, opt.MeanPowerN);
    sizeStats.meanRadius = momentRatioMean(radius, opt.MeanPowerM, opt.MeanPowerN);
    sizeStats.meanDiameter = momentRatioMean(diameter, opt.MeanPowerM, opt.MeanPowerN);
    sizeStats.areaWeightedMeanRadius = weightedMean(radius, area);
    sizeStats.areaWeightedMeanDiameter = weightedMean(diameter, area);
    sizeStats.hist = struct( ...
        'area', buildPdfHistogramData(area, round(opt.HistNumBins)), ...
        'radius', buildPdfHistogramData(radius, round(opt.HistNumBins)));
end

function value = maxOrNaN(x)
    if isempty(x)
        value = NaN;
    else
        value = max(x);
    end
end

function value = zeroIfEmpty(x)
    if isempty(x) || (isscalar(x) && isnan(x))
        value = 0;
    else
        value = x;
    end
end

function meanVal = weightedMean(x, w)
    x = x(:);
    w = w(:);
    valid = isfinite(x) & isfinite(w) & (w > 0);
    x = x(valid);
    w = w(valid);
    if isempty(x)
        meanVal = NaN;
        return;
    end
    meanVal = sum(x .* w) / sum(w);
end

function histData = buildPdfHistogramData(x, nBins)
    x = x(:);
    x = x(isfinite(x) & (x > 0));
    histData = struct('edges', [], 'centers', [], 'count', [], 'pdf', [], 'probability', []);
    if isempty(x)
        return;
    end

    nBins = max(1, round(nBins));
    xmin = min(x);
    xmax = max(x);
    if xmin == xmax
        pad = max(1e-12, abs(xmin) * 1e-6);
        edges = linspace(xmin - pad, xmax + pad, nBins + 1);
    else
        edges = linspace(xmin, xmax, nBins + 1);
    end
    count = histcounts(x, edges);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    widths = diff(edges);
    totalCount = sum(count);
    pdf = zeros(size(count));
    prob = zeros(size(count));
    if totalCount > 0
        prob = count ./ totalCount;
        pdf = count ./ (totalCount .* widths);
    end

    histData.edges = edges;
    histData.centers = centers;
    histData.count = count;
    histData.pdf = pdf;
    histData.probability = prob;
end

function net = buildPhaseSkeletonStats(phaseMask, phaseArea, dx, dy, boundary, enabled, toolbox)
    net = struct();
    net.enabled = enabled;
    net.available = toolbox.skeletonAvailable;
    net.skeletonLength = NaN;
    net.branchPoints = NaN;
    net.endPoints = NaN;
    net.branchDensity = NaN;
    net.note = '';

    if ~enabled
        net.note = 'Skeleton graph is disabled (EnableSkeletonGraph=false).';
        return;
    end
    if ~toolbox.skeletonAvailable
        net.note = 'Skeleton analysis requires Image Processing Toolbox (bwskel or bwmorph).';
        return;
    end
    if ~any(phaseMask(:))
        net.skeletonLength = 0;
        net.branchPoints = 0;
        net.endPoints = 0;
        net.branchDensity = 0;
        net.note = 'Phase is empty.';
        return;
    end

    [workMask, centerMask] = buildPeriodicMorphologyMask(phaseMask, boundary);
    skel = skeletonizeBinary(workMask, toolbox);
    branchMask = bwmorph(skel, 'branchpoints');
    endMask = bwmorph(skel, 'endpoints');

    net.skeletonLength = computeSkeletonLength(skel, centerMask, dx, dy);
    net.branchPoints = nnz(branchMask & centerMask);
    net.endPoints = nnz(endMask & centerMask);
    net.branchDensity = safeDivide(net.branchPoints, phaseArea);
    if ~all(centerMask(:))
        net.note = 'Periodic morphology metrics were computed on a tiled domain and cropped back to the center tile.';
    end
end

function thick = buildThicknessStats(matrixMask, dx, dy, boundary, toolbox)
    thick = struct('min', NaN, 'mean', NaN, 'p1', NaN, 'p5', NaN, ...
        'sampleCount', 0, 'note', '');
    if ~any(matrixMask(:))
        thick.note = 'Matrix phase is empty.';
        return;
    end
    if ~(toolbox.distanceAvailable && toolbox.skeletonAvailable)
        thick.note = 'Thickness analysis requires Image Processing Toolbox (bwdist plus bwskel/bwmorph).';
        return;
    end

    [workMask, centerMask] = buildPeriodicMorphologyMask(matrixMask, boundary);
    dist = bwdist(~workMask);
    scaleLen = effectivePixelLength(dx, dy);
    skel = skeletonizeBinary(workMask, toolbox);
    sampleMask = skel & centerMask;
    values = 2 .* dist(sampleMask) .* scaleLen;
    if isempty(values)
        values = 2 .* dist(workMask & centerMask) .* scaleLen;
        thick.note = 'Skeleton sampling was empty; thickness fell back to all matrix pixels in the center tile.';
    elseif anyPeriodicBoundary(boundary)
        thick.note = 'Thickness was sampled on the center tile of a periodic tiling.';
    end
    values = values(isfinite(values) & (values > 0));
    if isempty(values)
        return;
    end

    q = computeQuantiles(values, [0.01, 0.05]);
    thick.min = min(values);
    thick.mean = mean(values);
    thick.p1 = q(1);
    thick.p5 = q(2);
    thick.sampleCount = numel(values);
end

function frag = buildFragmentationStats(matrix, matrixSize, matrixTopology)
    frag = struct();
    frag.count = matrixTopology.beta0;
    frag.largestFraction = zeroIfNaN(matrix.largestComponentFraction);
    frag.sizeDist = matrixSize.hist.area;
end

function toolbox = detectImageToolboxAvailability()
    toolbox = struct();
    toolbox.bwskel = (exist('bwskel', 'file') == 2);
    toolbox.bwmorph = (exist('bwmorph', 'file') == 2);
    toolbox.distanceAvailable = (exist('bwdist', 'file') == 2);
    toolbox.skeletonAvailable = toolbox.bwskel || toolbox.bwmorph;
    toolbox.imageProcessingAvailable = toolbox.distanceAvailable && toolbox.skeletonAvailable;
end

function [workMask, centerMask] = buildPeriodicMorphologyMask(mask, boundary)
    [ny, nx] = size(mask);
    repY = 1;
    repX = 1;
    rowStart = 1;
    colStart = 1;
    if hasPeriodicY(boundary)
        repY = 3;
        rowStart = ny + 1;
    end
    if hasPeriodicX(boundary)
        repX = 3;
        colStart = nx + 1;
    end
    workMask = repmat(mask, repY, repX);
    centerMask = false(size(workMask));
    centerMask(rowStart:(rowStart + ny - 1), colStart:(colStart + nx - 1)) = true;
end

function skel = skeletonizeBinary(mask, toolbox)
    if toolbox.bwskel
        skel = bwskel(mask);
    else
        skel = bwmorph(mask, 'skel', Inf);
    end
    skel = logical(skel);
end

function len = computeSkeletonLength(skel, centerMask, dx, dy)
    [ny, nx] = size(skel);
    len = 0;
    [rows, cols] = find(skel & centerMask);
    for i = 1:numel(rows)
        r = rows(i);
        c = cols(i);
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue;
                end
                nr = r + dr;
                nc = c + dc;
                if nr < 1 || nr > ny || nc < 1 || nc > nx
                    continue;
                end
                if ~skel(nr, nc)
                    continue;
                end
                if centerMask(nr, nc)
                    if sub2ind([ny, nx], nr, nc) <= sub2ind([ny, nx], r, c)
                        continue;
                    end
                end
                len = len + neighborLength(dr, dc, dx, dy);
            end
        end
    end
end

function len = neighborLength(dr, dc, dx, dy)
    if dr == 0
        len = dy;
    elseif dc == 0
        len = dx;
    else
        len = hypot(dx, dy);
    end
end

function scaleLen = effectivePixelLength(dx, dy)
    scaleLen = sqrt(dx * dy);
end

function tf = anyPeriodicBoundary(boundary)
    tf = hasPeriodicX(boundary) || hasPeriodicY(boundary);
end

function out = buildEvolutionSeries(chunkFile, opt)
    axisName = resolveEvolutionAxis(opt);
    idx = loadChunkIndexCache(chunkFile);
    stepIndexAll = (1:numel(idx.timesteps)).';
    timestepAll = idx.timesteps(:);
    timeAll = nan(size(timestepAll));
    if strcmp(axisName, 'time')
        timeAll = mapTimestepsToPhysicalTime(timestepAll, opt, chunkFile);
        axisValueAll = timeAll;
    elseif ~isempty(opt.SlurmPath)
        timeAll = mapTimestepsToPhysicalTime(timestepAll, opt, chunkFile);
        axisValueAll = timestepAll;
    elseif strcmp(axisName, 'timestep')
        axisValueAll = timestepAll;
    else
        axisValueAll = stepIndexAll;
    end

    selectMask = true(size(axisValueAll));
    if ~isempty(opt.EvolutionRange)
        selectMask = selectMask & isfinite(axisValueAll) & ...
            (axisValueAll >= opt.EvolutionRange(1)) & (axisValueAll <= opt.EvolutionRange(2));
    end
    selected = find(selectMask);
    stride = max(1, round(opt.EvolutionStride));
    selected = selected(1:stride:end);

    out = emptyEvolutionStats();
    out.axis = axisName;
    if isempty(selected)
        return;
    end

    nSel = numel(selected);
    out.axisValue = axisValueAll(selected).';
    out.stepIndex = stepIndexAll(selected).';
    out.timestep = timestepAll(selected).';
    out.time = timeAll(selected).';

    out.geometry.phi = nan(1, nSel);
    out.geometry.interfaceLength = nan(1, nSel);
    out.geometry.specificInterface = nan(1, nSel);
    out.topology.pore.beta0 = nan(1, nSel);
    out.topology.pore.beta1 = nan(1, nSel);
    out.topology.pore.chi = nan(1, nSel);
    out.topology.matrix.beta0 = nan(1, nSel);
    out.topology.matrix.beta1 = nan(1, nSel);
    out.topology.matrix.chi = nan(1, nSel);
    out.connectivity.pore.largestFraction = nan(1, nSel);
    out.connectivity.pore.percolatesX = false(1, nSel);
    out.connectivity.pore.percolatesY = false(1, nSel);
    out.connectivity.matrix.largestFraction = nan(1, nSel);
    out.connectivity.matrix.percolatesX = false(1, nSel);
    out.connectivity.matrix.percolatesY = false(1, nSel);
    out.thickness.matrix.min = nan(1, nSel);
    out.thickness.matrix.mean = nan(1, nSel);
    out.thickness.matrix.p1 = nan(1, nSel);
    out.thickness.matrix.p5 = nan(1, nSel);
    out.fragmentation.matrix.count = nan(1, nSel);
    out.fragmentation.matrix.largestFraction = nan(1, nSel);

    for i = 1:nSel
        step = read_chunk_step_fast(chunkFile, 'SelectBy', 'Index', 'Index', out.stepIndex(i), ...
            'ProgressMode', 'off');
        snapStats = computeSnapshotPaperStats(step, chunkFile, opt);
        out.geometry.phi(i) = snapStats.geometry.phi;
        out.geometry.interfaceLength(i) = snapStats.geometry.interfaceLength;
        out.geometry.specificInterface(i) = snapStats.geometry.specificInterface;
        out.topology.pore.beta0(i) = snapStats.topology.pore.beta0;
        out.topology.pore.beta1(i) = snapStats.topology.pore.beta1;
        out.topology.pore.chi(i) = snapStats.topology.pore.chi;
        out.topology.matrix.beta0(i) = snapStats.topology.matrix.beta0;
        out.topology.matrix.beta1(i) = snapStats.topology.matrix.beta1;
        out.topology.matrix.chi(i) = snapStats.topology.matrix.chi;
        out.connectivity.pore.largestFraction(i) = snapStats.connectivity.pore.largestFraction;
        out.connectivity.pore.percolatesX(i) = snapStats.connectivity.pore.percolatesX;
        out.connectivity.pore.percolatesY(i) = snapStats.connectivity.pore.percolatesY;
        out.connectivity.matrix.largestFraction(i) = snapStats.connectivity.matrix.largestFraction;
        out.connectivity.matrix.percolatesX(i) = snapStats.connectivity.matrix.percolatesX;
        out.connectivity.matrix.percolatesY(i) = snapStats.connectivity.matrix.percolatesY;
        out.thickness.matrix.min(i) = snapStats.thickness.matrix.min;
        out.thickness.matrix.mean(i) = snapStats.thickness.matrix.mean;
        out.thickness.matrix.p1(i) = snapStats.thickness.matrix.p1;
        out.thickness.matrix.p5(i) = snapStats.thickness.matrix.p5;
        out.fragmentation.matrix.count(i) = snapStats.fragmentation.matrix.count;
        out.fragmentation.matrix.largestFraction(i) = snapStats.fragmentation.matrix.largestFraction;
    end

    out.transition.deltaPhi = diffWithLeadingNaN(out.geometry.phi);
    out.transition.deltaSpecificInterface = diffWithLeadingNaN(out.geometry.specificInterface);
    out.transition.pore.deltaBeta0 = diffWithLeadingNaN(out.topology.pore.beta0);
    out.transition.pore.deltaBeta1 = diffWithLeadingNaN(out.topology.pore.beta1);
    out.transition.pore.deltaLargestFraction = diffWithLeadingNaN(out.connectivity.pore.largestFraction);
    out.transition.matrix.deltaBeta0 = diffWithLeadingNaN(out.topology.matrix.beta0);
    out.transition.matrix.deltaBeta1 = diffWithLeadingNaN(out.topology.matrix.beta1);
    out.transition.matrix.deltaLargestFraction = diffWithLeadingNaN(out.connectivity.matrix.largestFraction);
end

function stats = computeSnapshotPaperStats(step, chunkFile, opt)
    data = step.data;
    col = step.colIndex;
    [x, y, ncount] = extractRequiredColumns(data, col, opt.NcountVar);
    [x, y] = scaleCoordinates(x, y, opt.CoordScale);
    [ncountGrid, validMask, xCenters, yCenters] = buildGrid(x, y, ncount);
    [dx, dy] = resolveGridSpacing(opt, chunkFile, xCenters, yCenters);
    poreMask = validMask & (ncountGrid < opt.ThresholdN);
    matrixMask = validMask & ~poreMask;
    pore = analyzePhase('pore', poreMask, validMask, xCenters, yCenters, dx, dy, opt);
    matrix = analyzePhase('matrix', matrixMask, validMask, xCenters, yCenters, dx, dy, opt);

    cellArea = dx * dy;
    validArea = nnz(validMask) * cellArea;
    poreArea = nnz(poreMask) * cellArea;
    matrixArea = nnz(matrixMask) * cellArea;
    globalStats = struct();
    globalStats.validArea = validArea;
    globalStats.poreArea = poreArea;
    globalStats.matrixArea = matrixArea;
    globalStats.porosity = safeDivide(poreArea, validArea);
    globalStats.matrixFraction = safeDivide(matrixArea, validArea);
    globalStats.interfaceLength = pore.totalInterfacePerimeter;
    globalStats.specificInterfaceBulk = safeDivide(globalStats.interfaceLength, validArea);
    globalStats.specificInterfacePore = safeDivide(globalStats.interfaceLength, poreArea);
    globalStats.specificInterfaceMatrix = safeDivide(globalStats.interfaceLength, matrixArea);
    stats = buildPaperStats(globalStats, pore, matrix, poreMask, matrixMask, validMask, ...
        xCenters, yCenters, dx, dy, step.timestep, chunkFile, opt);
end

function axisName = resolveEvolutionAxis(opt)
    axisName = opt.EvolutionSelectBy;
    if isempty(axisName)
        axisName = lower(strtrim(opt.SelectBy));
    end
    if strcmp(axisName, 'timestep')
        return;
    end
    if strcmp(axisName, 'time')
        return;
    end
    axisName = 'index';
end

function idx = loadChunkIndexCache(chunkFile)
    cachePath = [chunkFile, '.stepidx.mat'];
    if ~exist(cachePath, 'file')
        read_chunk_step_fast(chunkFile, 'SelectBy', 'Index', 'Index', 1, 'ProgressMode', 'off');
    end
    S = load(cachePath, 'idx');
    idx = S.idx;
end

function timeAll = mapTimestepsToPhysicalTime(timestepAll, opt, chunkFile)
    slurmPath = opt.SlurmPath;
    if isempty(slurmPath)
        slurmPath = get_slurm_txt_fullpath(fileparts(chunkFile));
    end
    if isempty(slurmPath) || ~exist(slurmPath, 'file')
        error('analyze_chunk_network2d:MissingEvolutionSlurmPath', ...
            'EvolutionSelectBy=Time requires a valid SlurmPath.');
    end

    S = read_slurm_stepcpu(slurmPath);
    m = round(opt.SlurmModuleIndex);
    if m < 1 || m > numel(S.modules)
        error('analyze_chunk_network2d:BadEvolutionSlurmModule', ...
            'SlurmModuleIndex out of range for evolution time mapping.');
    end
    scol = S.modules(m).colIndex;
    if ~isfield(scol, 'Time') || ~isfield(scol, 'Step')
        error('analyze_chunk_network2d:BadEvolutionSlurmColumns', ...
            'Slurm block must contain Time and Step columns for evolution mapping.');
    end

    slurmTimeAll = S.modules(m).data(:, scol.Time);
    slurmStepAll = S.modules(m).data(:, scol.Step);
    valid = isfinite(slurmTimeAll) & isfinite(slurmStepAll);
    slurmTime = slurmTimeAll(valid);
    slurmStep = slurmStepAll(valid);
    timeAll = nan(size(timestepAll));
    for i = 1:numel(timestepAll)
        [~, idxMin] = min(abs(slurmStep - timestepAll(i)));
        timeAll(i) = slurmTime(idxMin);
    end
end

function d = diffWithLeadingNaN(x)
    x = x(:).';
    d = nan(size(x));
    if numel(x) >= 2
        d(2:end) = diff(x);
    end
end

function out = emptyEvolutionStats()
    out = struct();
    out.axis = '';
    out.axisValue = [];
    out.stepIndex = [];
    out.timestep = [];
    out.time = [];
    out.geometry = struct('phi', [], 'interfaceLength', [], 'specificInterface', []);
    out.topology = struct( ...
        'pore', struct('beta0', [], 'beta1', [], 'chi', []), ...
        'matrix', struct('beta0', [], 'beta1', [], 'chi', []));
    out.connectivity = struct( ...
        'pore', struct('largestFraction', [], 'percolatesX', [], 'percolatesY', []), ...
        'matrix', struct('largestFraction', [], 'percolatesX', [], 'percolatesY', []));
    out.thickness = struct('matrix', struct('min', [], 'mean', [], 'p1', [], 'p5', []));
    out.fragmentation = struct('matrix', struct('count', [], 'largestFraction', []));
    out.transition = struct( ...
        'deltaPhi', [], ...
        'deltaSpecificInterface', [], ...
        'pore', struct('deltaBeta0', [], 'deltaBeta1', [], 'deltaLargestFraction', []), ...
        'matrix', struct('deltaBeta0', [], 'deltaBeta1', [], 'deltaLargestFraction', []));
end

function fig = plotPhaseGrid(xCenters, yCenters, phaseGrid, timestep, thresholdN, xRange, yRange)
    fig = figure('Color', 'w', 'Name', '2D Network Phase Map');
    ax = axes('Parent', fig);
    imagesc(ax, xCenters, yCenters, phaseGrid);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    applyPlotRanges(ax, xCenters, yCenters, xRange, yRange);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.12, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
    colormap(ax, [0.75 0.75 0.75; 0.20 0.45 0.80; 0.88 0.52 0.22]);
    setAxesCLim(ax, [-1 1]);
    cb = colorbar(ax, 'Ticks', [-1, 0, 1], 'TickLabels', {'Invalid', 'Matrix', 'Pore'});
    ylabel(cb, 'Phase');
    xlabel(ax, 'x', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'y', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, sprintf('Phase map (Ncount < %g is pore) @ timestep %g', thresholdN, timestep), ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = plotLabelGrid(xCenters, yCenters, labelGrid, validMask, ttl, xRange, yRange)
    fig = figure('Color', 'w', 'Name', ttl);
    ax = axes('Parent', fig);
    z = labelGrid;
    z(~validMask) = NaN;
    imagesc(ax, xCenters, yCenters, z);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    applyPlotRanges(ax, xCenters, yCenters, xRange, yRange);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.12, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
    cmap = lines(max(max(labelGrid(:)), 1) + 1);
    colormap(ax, cmap);
    cb = colorbar(ax);
    ylabel(cb, 'Component Label');
    xlabel(ax, 'x', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'y', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, ttl, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = plotConnectivityHighlights(xCenters, yCenters, validMask, pore, matrix, timestep, xRange, yRange)
    fig = figure('Color', 'w', 'Name', 'Directional Connectivity Highlights');

    subplot(2, 2, 1);
    drawConnectivityPanel(xCenters, yCenters, validMask, pore.mask, pore.labelGrid, ...
        pore.connectivity.x.componentIds, 'Pore connected in x', xRange, yRange);

    subplot(2, 2, 2);
    drawConnectivityPanel(xCenters, yCenters, validMask, pore.mask, pore.labelGrid, ...
        pore.connectivity.y.componentIds, 'Pore connected in y', xRange, yRange);

    subplot(2, 2, 3);
    drawConnectivityPanel(xCenters, yCenters, validMask, matrix.mask, matrix.labelGrid, ...
        matrix.connectivity.x.componentIds, 'Matrix connected in x', xRange, yRange);

    subplot(2, 2, 4);
    drawConnectivityPanel(xCenters, yCenters, validMask, matrix.mask, matrix.labelGrid, ...
        matrix.connectivity.y.componentIds, 'Matrix connected in y', xRange, yRange);

    annotation(fig, 'textbox', [0.28 0.95 0.44 0.04], 'String', ...
        sprintf('Directional connectivity highlights @ timestep %g', timestep), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
end

function drawConnectivityPanel(xCenters, yCenters, validMask, phaseMask, labelGrid, componentIds, ttl, xRange, yRange)
    ax = gca;
    highlightMask = false(size(labelGrid));
    if ~isempty(componentIds)
        highlightMask = ismember(labelGrid, componentIds);
    end

    z = nan(size(labelGrid));
    z(validMask & ~phaseMask) = 0;
    z(validMask & phaseMask & ~highlightMask) = 1;
    z(validMask & highlightMask) = 2;

    imagesc(ax, xCenters, yCenters, z);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    applyPlotRanges(ax, xCenters, yCenters, xRange, yRange);
    box(ax, 'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.10, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 11);
    colormap(ax, [1.00 1.00 1.00; 0.77 0.84 0.93; 0.88 0.42 0.05]);
    setAxesCLim(ax, [0 2]);
    xlabel(ax, 'x', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax, 'y', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax, ttl, 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    if ~any(highlightMask(:))
        xLimNow = get(ax, 'XLim');
        yLimNow = get(ax, 'YLim');
        text(ax, mean(xLimNow), mean(yLimNow), 'No connected component', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', [0.25 0.25 0.25], 'FontName', 'Times New Roman', ...
            'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1]);
    end
end

function applyPlotRanges(ax, xCenters, yCenters, xRange, yRange)
    if ~isempty(xRange)
        xlim(ax, normalizePlotRange(xRange, xCenters));
    end
    if ~isempty(yRange)
        ylim(ax, normalizePlotRange(yRange, yCenters));
    end
end

function rangeOut = normalizePlotRange(rangeIn, centers)
    rangeOut = double(rangeIn(:).');
    if isempty(rangeOut) || (rangeOut(2) > rangeOut(1))
        return;
    end
    padding = inferPlotPadding(centers);
    rangeOut = [rangeOut(1) - padding, rangeOut(2) + padding];
end

function padding = inferPlotPadding(centers)
    centers = double(centers(:));
    centers = centers(isfinite(centers));
    if numel(centers) >= 2
        d = diff(sort(centers));
        d = d(d > 0);
        if ~isempty(d)
            padding = 0.5 * min(d);
            return;
        end
    end
    padding = 0.5;
end

function fig = plotCountDistribution(dist, phaseName)
    fig = [];
    if isempty(dist.count)
        return;
    end
    fig = figure('Color', 'w', 'Name', [phaseName, ' Count Distribution']);
    ax = axes('Parent', fig);
    bar(ax, dist.centers, dist.count, 1.0, 'FaceColor', [0.35 0.6 0.85], 'EdgeColor', 'none');
    xlabel(ax, 'Equivalent Diameter', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'Count', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, [phaseName, ' size count distribution'], ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.16, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
end

function fig = plotPositionDistribution(dist, phaseName, axisName)
    fig = [];
    if isempty(dist.meanEquivDiameter)
        return;
    end
    valid = (dist.count > 0) & isfinite(dist.meanEquivDiameter);
    if ~any(valid)
        return;
    end
    fig = figure('Color', 'w', 'Name', [phaseName, ' Position Distribution ', upper(axisName)]);
    ax = axes('Parent', fig);
    plot(ax, dist.centers(valid), dist.meanEquivDiameter(valid), 'o-', ...
        'LineWidth', 1.5, 'Color', [0.12 0.40 0.75], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.12 0.40 0.75]);
    xlabel(ax, axisName, 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'Mean Equivalent Diameter', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, sprintf('%s mean size vs %s bin', phaseName, axisName), ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.16, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
end

function fig = plotDirectionalProfile(profile, timestep)
    fig = figure('Color', 'w', 'Name', ['Directional Profile ', upper(profile.axis)]);

    ax1 = subplot(3, 1, 1, 'Parent', fig);
    plot(ax1, profile.centers, profile.porosity, 'o-', ...
        'LineWidth', 1.5, 'Color', [0.10 0.42 0.78], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.10 0.42 0.78]);
    ylabel(ax1, 'Porosity', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax1, sprintf('Pore-structure profile along %s @ timestep %g', profile.axis, timestep), ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    styleProfileAxis(ax1);

    ax2 = subplot(3, 1, 2, 'Parent', fig);
    plot(ax2, profile.centers, profile.specificInterfaceBulk, 's-', ...
        'LineWidth', 1.5, 'Color', [0.85 0.38 0.08], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.85 0.38 0.08]);
    ylabel(ax2, 'Specific Interface', 'FontName', 'Times New Roman', 'FontSize', 12);
    styleProfileAxis(ax2);

    ax3 = subplot(3, 1, 3, 'Parent', fig);
    yyaxis(ax3, 'left');
    plot(ax3, profile.centers, profile.connectivityFraction, '^-', ...
        'LineWidth', 1.5, 'Color', [0.15 0.55 0.20], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.15 0.55 0.20]);
    ylabel(ax3, 'Conn. Fraction', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylim(ax3, [0 1]);
    yyaxis(ax3, 'right');
    stairs(ax3, profile.centers, double(profile.connectivityFlag), '--', ...
        'LineWidth', 1.4, 'Color', [0.55 0.10 0.10]);
    ylabel(ax3, 'Connected (0/1)', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylim(ax3, [0 1]);
    xlabel(ax3, profile.axis, 'FontName', 'Times New Roman', 'FontSize', 13);
    legend(ax3, {'Connected area fraction', 'Perp. connected flag'}, 'Location', 'best');
    styleProfileAxis(ax3);
end

function fig = plotEvolutionSeries(evolution)
    fig = [];
    if isempty(evolution.axisValue)
        return;
    end

    x = evolution.axisValue;
    if any(~isfinite(x))
        x = evolution.stepIndex;
        xLabel = 'index';
    else
        xLabel = evolution.axis;
    end

    fig = figure('Color', 'w', 'Name', 'Network2d Evolution');

    ax1 = subplot(2, 2, 1, 'Parent', fig);
    plot(ax1, x, evolution.geometry.phi, 'o-', 'LineWidth', 1.5, ...
        'Color', [0.10 0.42 0.78], 'MarkerFaceColor', [0.10 0.42 0.78]);
    xlabel(ax1, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax1, '\phi', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax1, 'Porosity Evolution', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    styleProfileAxis(ax1);

    ax2 = subplot(2, 2, 2, 'Parent', fig);
    plot(ax2, x, evolution.geometry.specificInterface, 's-', 'LineWidth', 1.5, ...
        'Color', [0.85 0.38 0.08], 'MarkerFaceColor', [0.85 0.38 0.08]);
    xlabel(ax2, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax2, 'L / area', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax2, 'Specific Interface Evolution', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    styleProfileAxis(ax2);

    ax3 = subplot(2, 2, 3, 'Parent', fig);
    plot(ax3, x, evolution.connectivity.pore.largestFraction, 'o-', 'LineWidth', 1.5, ...
        'Color', [0.15 0.55 0.20], 'MarkerFaceColor', [0.15 0.55 0.20]);
    hold(ax3, 'on');
    plot(ax3, x, evolution.connectivity.matrix.largestFraction, 'd-', 'LineWidth', 1.5, ...
        'Color', [0.55 0.10 0.10], 'MarkerFaceColor', [0.55 0.10 0.10]);
    xlabel(ax3, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax3, 'Largest Fraction', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax3, 'Largest Cluster Fraction', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    legend(ax3, {'Pore', 'Matrix'}, 'Location', 'best');
    styleProfileAxis(ax3);

    ax4 = subplot(2, 2, 4, 'Parent', fig);
    plot(ax4, x, evolution.topology.pore.beta0, 'o-', 'LineWidth', 1.3, ...
        'Color', [0.10 0.42 0.78], 'MarkerFaceColor', [0.10 0.42 0.78]);
    hold(ax4, 'on');
    plot(ax4, x, evolution.topology.pore.beta1, 's-', 'LineWidth', 1.3, ...
        'Color', [0.85 0.38 0.08], 'MarkerFaceColor', [0.85 0.38 0.08]);
    plot(ax4, x, evolution.topology.matrix.beta0, 'd--', 'LineWidth', 1.3, ...
        'Color', [0.15 0.55 0.20], 'MarkerFaceColor', [0.15 0.55 0.20]);
    plot(ax4, x, evolution.topology.matrix.beta1, '^-', 'LineWidth', 1.3, ...
        'Color', [0.55 0.10 0.10], 'MarkerFaceColor', [0.55 0.10 0.10]);
    xlabel(ax4, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax4, 'Betti Number', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax4, 'Topology Evolution', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    legend(ax4, {'pore \beta_0', 'pore \beta_1', 'matrix \beta_0', 'matrix \beta_1'}, 'Location', 'best');
    styleProfileAxis(ax4);
end

function styleProfileAxis(ax)
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.16, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
end

function plots = emptyPlotsStruct()
    plots = struct( ...
        'phaseFig', [], ...
        'poreLabelFig', [], ...
        'matrixLabelFig', [], ...
        'connectivityFig', [], ...
        'evolutionFig', [], ...
        'profileFigX', [], ...
        'profileFigY', [], ...
        'poreCountFig', [], ...
        'poreProbFig', [], ...
        'poreCdfFig', [], ...
        'matrixCountFig', [], ...
        'matrixProbFig', [], ...
        'matrixCdfFig', [], ...
        'porePositionFigX', [], ...
        'porePositionFigY', [], ...
        'matrixPositionFigX', [], ...
        'matrixPositionFigY', []);
end

function comp = emptyComponentStats()
    comp = struct( ...
        'label', zeros(0, 1), ...
        'cellCount', zeros(0, 1), ...
        'area', zeros(0, 1), ...
        'equivDiameter', zeros(0, 1), ...
        'centroidX', zeros(0, 1), ...
        'centroidY', zeros(0, 1), ...
        'bboxWidth', zeros(0, 1), ...
        'bboxHeight', zeros(0, 1), ...
        'perimeterOpen', zeros(0, 1), ...
        'interfacePerimeter', zeros(0, 1), ...
        'shapeFactor', zeros(0, 1), ...
        'touchesLeft', false(0, 1), ...
        'touchesRight', false(0, 1), ...
        'touchesBottom', false(0, 1), ...
        'touchesTop', false(0, 1), ...
        'wrapX', false(0, 1), ...
        'wrapY', false(0, 1), ...
        'percolatesX', false(0, 1), ...
        'percolatesY', false(0, 1));
end

function [nr, nc, hasNeighbor, seamX, seamY] = neighborAt(r, c, dirId, ny, nx, boundary)
    [nr, nc, hasNeighbor, shiftX, shiftY] = neighborAtWithShift(r, c, dirId, ny, nx, boundary);
    seamX = (shiftX ~= 0);
    seamY = (shiftY ~= 0);
end

function [nr, nc, hasNeighbor, shiftX, shiftY] = neighborAtWithShift(r, c, dirId, ny, nx, boundary)
    nr = r;
    nc = c;
    hasNeighbor = true;
    shiftX = 0;
    shiftY = 0;

    switch dirId
        case 1
            nc = c - 1;
            if nc < 1
                if hasPeriodicX(boundary)
                    nc = nx;
                    shiftX = -1;
                else
                    hasNeighbor = false;
                end
            end
        case 2
            nc = c + 1;
            if nc > nx
                if hasPeriodicX(boundary)
                    nc = 1;
                    shiftX = 1;
                else
                    hasNeighbor = false;
                end
            end
        case 3
            nr = r - 1;
            if nr < 1
                if hasPeriodicY(boundary)
                    nr = ny;
                    shiftY = -1;
                else
                    hasNeighbor = false;
                end
            end
        case 4
            nr = r + 1;
            if nr > ny
                if hasPeriodicY(boundary)
                    nr = 1;
                    shiftY = 1;
                else
                    hasNeighbor = false;
                end
            end
        otherwise
            hasNeighbor = false;
    end
end

function len = sideLengthByDir(dirId, dx, dy)
    if dirId == 1 || dirId == 2
        len = dy;
    else
        len = dx;
    end
end

function tf = hasPeriodicX(boundary)
    tf = strcmp(boundary, 'periodic-x') || strcmp(boundary, 'periodic-xy');
end

function tf = hasPeriodicY(boundary)
    tf = strcmp(boundary, 'periodic-y') || strcmp(boundary, 'periodic-xy');
end

function info = buildSelectionInfo(step, opt)
    mode = lower(strtrim(opt.SelectBy));
    info = struct();
    info.mode = mode;
    info.stepIndex = step.stepIndex;
    info.mappedTimeStep = step.timestep;

    switch mode
        case 'index'
            info.requestedIndex = opt.Index;
        case 'timestep'
            info.requestedTimeStep = opt.TimeStep;
        case 'time'
            info.requestedTime = opt.Time;
        otherwise
            % read_chunk_step_fast has already validated selector.
    end
end

function y = safeDivide(a, b)
    if isempty(b) || ~isfinite(b) || b == 0
        y = NaN;
    else
        y = a ./ b;
    end
end

function y = safeDivideOrZero(a, b)
    if isempty(b) || ~isfinite(b) || b == 0
        if isempty(a) || ~isfinite(a)
            y = NaN;
        elseif a == 0
            y = 0;
        else
            y = NaN;
        end
    else
        y = a ./ b;
    end
end

function y = zeroIfNaN(x)
    if isempty(x) || (isscalar(x) && isnan(x))
        y = 0;
    else
        y = x;
    end
end

function setAxesCLim(ax, limits)
    if nargin < 2 || isempty(limits)
        return;
    end
    set(ax, 'CLim', limits);
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
