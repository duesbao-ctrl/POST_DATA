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
    p.addParameter('DiameterRange', [], @isnumeric);
    p.addParameter('DiameterHistBinSize', [], @isnumeric);
    p.addParameter('DiameterEmptyBinMode', 'zero', @isTextScalar);
    p.addParameter('DiameterPlotStyle', 'bar', @isTextScalar);
    p.addParameter('DiameterFitTypes', {'powerlaw', 'gamma', 'lognormal'});
    p.addParameter('HistScale', 'linear', @isTextScalar);
    p.addParameter('MeanPowerM', 1, @isnumeric);
    p.addParameter('MeanPowerN', 0, @isnumeric);
    p.addParameter('GeometryMode', 'cutcell', @isTextScalar);
    p.addParameter('CutCellMethod', 'plic', @isTextScalar);
    p.addParameter('CutCellFallback', 'binary', @isTextScalar);
    p.addParameter('CutCellPlotRefinement', 4, @isnumeric);
    p.addParameter('EnableSkeletonGraph', false, @islogical);
    p.addParameter('EnableEvolution', false, @islogical);
    p.addParameter('EvolutionSelectBy', '', @isTextScalar);
    p.addParameter('EvolutionRange', [], @isnumeric);
    p.addParameter('EvolutionStride', 1, @isnumeric);
    p.addParameter('MakePlots', true, @islogical);
    p.addParameter('CancelCallback', @() false, @(x) isa(x, 'function_handle'));
    p.addParameter('ProgressCallback', @(fraction, message) [], @(x) isa(x, 'function_handle'));
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
    opt.DiameterEmptyBinMode = pd_distribution_normalize_empty_bin_mode( ...
        opt.DiameterEmptyBinMode);
    opt.DiameterPlotStyle = lower(strtrim(toChar(opt.DiameterPlotStyle)));
    opt.DiameterFitTypes = pd_distribution_normalize_fit_types( ...
        opt.DiameterFitTypes, 'DiameterFitTypes', ...
        'analyze_chunk_network2d:BadFitTypes');
    opt.HistScale = lower(strtrim(toChar(opt.HistScale)));
    opt.GeometryMode = normalizeGeometryMode(opt.GeometryMode);
    opt.CutCellMethod = lower(strtrim(toChar(opt.CutCellMethod)));
    opt.CutCellFallback = lower(strtrim(toChar(opt.CutCellFallback)));

    validateInputs(opt);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex, ...
                    'ProgressMode', opt.ProgressMode, ...
                    'CancelCallback', opt.CancelCallback, ...
                    'ProgressCallback', opt.ProgressCallback};

    step = read_chunk_step_fast(chunkFile, selectorArgs{:});
    selectionInfo = buildSelectionInfo(step, opt);
    data = step.data;
    col = step.colIndex;

    [x, y, ncount] = pd_network_extract_columns(data, col, opt.NcountVar);
    [x, y] = pd_network_scale_coordinates(x, y, opt.CoordScale);
    [ncountGrid, validMask, xCenters, yCenters, presentCount] = ...
        pd_network_build_grid(x, y, ncount);
    [dx, dy, spacingSource] = pd_network_resolve_spacing(opt, chunkFile, xCenters, yCenters);

    [poreMask, matrixMask] = pd_network_classify_phases( ...
        ncountGrid, validMask, opt.ThresholdN);
    cutCell = pd_network_build_geometry(ncountGrid, validMask, poreMask, matrixMask, ...
        xCenters, yCenters, dx, dy, opt);
    phaseGrid = cutCell.phaseGrid;

    pore = analyzePhase('pore', poreMask, validMask, xCenters, yCenters, dx, dy, opt, cutCell.pore);
    matrix = analyzePhase('matrix', matrixMask, validMask, xCenters, yCenters, dx, dy, opt, cutCell.matrix);
    profile = pd_network_build_directional_profiles( ...
        poreMask, validMask, xCenters, yCenters, dx, dy, opt, pore, cutCell);

    cellArea = dx * dy;
    numValidCells = nnz(validMask);
    validArea = numValidCells * cellArea;
    poreCells = nnz(poreMask);
    matrixCells = nnz(matrixMask);
    poreArea = pore.area;
    matrixArea = matrix.area;
    interfaceLength = cutCell.interfaceLength;

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

    stats = pd_network_build_advanced_stats( ...
        globalStats, pore, matrix, poreMask, matrixMask, validMask, ...
        xCenters, yCenters, dx, dy, step.timestep, chunkFile, opt, cutCell);

    plots = pd_render_network2d_snapshot(xCenters, yCenters, phaseGrid, ...
        validMask, pore, matrix, profile, cutCell, step.timestep, opt);

    if opt.EnableEvolution
        stats.evolution = pd_network_build_evolution(chunkFile, opt, ...
            @(step) computeSnapshotPaperStats(step, chunkFile, opt));
        if opt.MakePlots
            plots.evolutionFig = pd_render_network2d_evolution(stats.evolution);
        end
    else
        stats.evolution = pd_network_empty_evolution();
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
    out.diameterRange = opt.DiameterRange;
    out.diameterEmptyBinMode = opt.DiameterEmptyBinMode;
    out.diameterPlotStyle = opt.DiameterPlotStyle;
    out.geometryMode = opt.GeometryMode;
    out.cutCell = cutCell;
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
    validateHistScale(opt.HistScale, 'analyze_chunk_network2d:BadHistScale');
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
    validatePositiveScalarOrEmpty(opt.DiameterHistBinSize, 'DiameterHistBinSize');
    validateRangeOrEmpty(opt.DiameterRange, 'DiameterRange');
    validateEmptyBinMode(opt.DiameterEmptyBinMode);
    validateDiameterPlotStyle(opt.DiameterPlotStyle);
    validateRangeOrEmpty(opt.PositionRangeX, 'PositionRangeX');
    validateRangeOrEmpty(opt.PositionRangeY, 'PositionRangeY');
    validateRangeOrEmpty(opt.PlotRangeX, 'PlotRangeX');
    validateRangeOrEmpty(opt.PlotRangeY, 'PlotRangeY');
    validateRangeOrEmpty(opt.ProfileRangeX, 'ProfileRangeX');
    validateRangeOrEmpty(opt.ProfileRangeY, 'ProfileRangeY');
    validateRangeOrEmpty(opt.EvolutionRange, 'EvolutionRange');
    validateCutCellOptions(opt);
end

function validateCutCellOptions(opt)
    validMode = {'cutcell', 'original'};
    if ~any(strcmp(opt.GeometryMode, validMode))
        error('analyze_chunk_network2d:BadGeometryMode', ...
            'GeometryMode must be cutcell or original.');
    end
    validMethod = {'plic'};
    if ~any(strcmp(opt.CutCellMethod, validMethod))
        error('analyze_chunk_network2d:BadCutCellMethod', ...
            'CutCellMethod must be plic.');
    end
    validFallback = {'binary', 'error'};
    if ~any(strcmp(opt.CutCellFallback, validFallback))
        error('analyze_chunk_network2d:BadCutCellFallback', ...
            'CutCellFallback must be binary or error.');
    end
    if ~(isscalar(opt.CutCellPlotRefinement) && isnumeric(opt.CutCellPlotRefinement) && ...
            isfinite(opt.CutCellPlotRefinement) && opt.CutCellPlotRefinement >= 1)
        error('analyze_chunk_network2d:BadCutCellPlotRefinement', ...
            'CutCellPlotRefinement must be a positive scalar.');
    end
end

function mode = normalizeGeometryMode(value)
    mode = lower(strtrim(toChar(value)));
    switch mode
        case {'cutcell', 'cut-cell', 'cut_cell', 'plic'}
            mode = 'cutcell';
        case {'original', 'binary', 'raw', 'grid', 'cell'}
            mode = 'original';
    end
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

function validateEmptyBinMode(mode)
    valid = {'zero', 'nan', 'remove'};
    if ~any(strcmp(mode, valid))
        error('analyze_chunk_network2d:BadEmptyBinMode', ...
            'DiameterEmptyBinMode must be zero/nan/remove.');
    end
end

function validateDiameterPlotStyle(plotStyle)
    valid = {'bar', 'scatter'};
    if ~any(strcmp(plotStyle, valid))
        error('analyze_chunk_network2d:BadDiameterPlotStyle', ...
            'DiameterPlotStyle must be bar or scatter.');
    end
end

function validateHistScale(histScale, errorId)
    valid = {'linear', 'semilogx', 'semilogy', 'loglog', 'semilog', 'log'};
    if ~any(strcmpi(strtrim(histScale), valid))
        error(errorId, 'HistScale must be linear/semilogx/semilogy/loglog.');
    end
end

function phase = analyzePhase(name, phaseMask, validMask, xCenters, yCenters, dx, dy, opt, phaseGeom)
    if nargin < 9
        phaseGeom = [];
    end
    [labelGrid, wrapXComp, wrapYComp] = ...
        pd_network_label_components(phaseMask, validMask, opt.Boundary);
    numComponents = max(labelGrid(:));
    [comp, sizeRank, componentStatus] = pd_network_compute_component_stats( ...
        labelGrid, phaseMask, validMask, xCenters, yCenters, dx, dy, ...
        opt.Boundary, wrapXComp, wrapYComp, phaseGeom);

    if isstruct(phaseGeom) && isfield(phaseGeom, 'areaGrid')
        area = sum(phaseGeom.areaGrid(:));
    else
        area = nnz(phaseMask) * dx * dy;
    end
    validArea = nnz(validMask) * dx * dy;
    largestArea = 0;
    if ~isempty(comp.area)
        largestArea = max(comp.area);
    end

    hasMissingCells = any(~validMask(:));
    topology = pd_network_compute_topology(phaseMask, opt.Boundary, numComponents, hasMissingCells);
    holeCount = topology.beta1;
    eulerCharacteristic = topology.chi;
    topologyNote = topology.note;

    [diamForSize, areaForSize, sizeMask] = pd_network_filter_component_sizes( ...
        comp.equivDiameter, comp.area, opt.DiameterRange);
    areaStats = pd_stats_summary(areaForSize);
    diamStats = pd_stats_summary(diamForSize, opt.MeanPowerM, opt.MeanPowerN);
    equivDiameterDistribution = buildHistogramDistribution(diamForSize, opt.DiameterRange, ...
        opt.DiameterHistBinSize, opt.DiameterFitTypes, opt.DiameterEmptyBinMode);
    positionDistribution = struct();
    positionDistribution.x = buildPositionDistribution(comp.centroidX(sizeMask), diamForSize, ...
        round(opt.PositionNumBins), opt.PositionRangeX, dx, opt.MeanPowerM, opt.MeanPowerN);
    positionDistribution.y = buildPositionDistribution(comp.centroidY(sizeMask), diamForSize, ...
        round(opt.PositionNumBins), opt.PositionRangeY, dy, opt.MeanPowerM, opt.MeanPowerN);
    connectivity = pd_network_build_directional_connectivity(comp, opt.Boundary, area);

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
    phase.percolatesX = componentStatus.percolatesX;
    phase.percolatesY = componentStatus.percolatesY;
    phase.wrapsX = componentStatus.wrapsX;
    phase.wrapsY = componentStatus.wrapsY;
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
    phase.sizeRank = sizeRank;
    phase.areaStats = areaStats;
    phase.equivDiameterStats = diamStats;
    phase.equivDiameterDistribution = equivDiameterDistribution;
    phase.positionDistribution = positionDistribution;
    phase.connectivity = connectivity;
end

function dist = buildHistogramDistribution(x, diameterRange, binSize, fitTypes, emptyBinMode)
    shared = pd_distribution_build(x, diameterRange, binSize, fitTypes, emptyBinMode);
    dist = struct('edges', shared.edges, 'centers', shared.centers, ...
        'count', shared.count, 'probability', shared.probability, ...
        'rawCenters', shared.rawCenters, 'rawCount', shared.rawCount, ...
        'rawProbability', shared.rawProbability, 'keptBins', shared.keptBins, ...
        'emptyBinMode', shared.emptyBinMode, 'cdfX', shared.cdfX, ...
        'cdfCount', shared.cdfCount, 'cdfProbability', shared.cdfProbability, ...
        'binSize', shared.binSize, 'fit', shared.fit);
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
        binId(i) = pd_network_locate_bin(position(i), edges);
    end
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    meanSize = nan(size(centers));
    count = zeros(size(centers));
    for i = 1:numel(centers)
        idx = (binId == i);
        count(i) = sum(idx);
        if count(i) > 0
            meanSize(i) = pd_stats_moment_ratio(sizeValue(idx), meanPowerM, meanPowerN);
        end
    end

    out.edges = edges;
    out.centers = centers;
    out.meanEquivDiameter = meanSize;
    out.count = count;
end

function stats = computeSnapshotPaperStats(step, chunkFile, opt)
    data = step.data;
    col = step.colIndex;
    [x, y, ncount] = pd_network_extract_columns(data, col, opt.NcountVar);
    [x, y] = pd_network_scale_coordinates(x, y, opt.CoordScale);
    [ncountGrid, validMask, xCenters, yCenters] = ...
        pd_network_build_grid(x, y, ncount);
    [dx, dy] = pd_network_resolve_spacing(opt, chunkFile, xCenters, yCenters);
    [poreMask, matrixMask] = pd_network_classify_phases( ...
        ncountGrid, validMask, opt.ThresholdN);
    cutCell = pd_network_build_geometry(ncountGrid, validMask, poreMask, matrixMask, ...
        xCenters, yCenters, dx, dy, opt);
    pore = analyzePhase('pore', poreMask, validMask, xCenters, yCenters, dx, dy, opt, cutCell.pore);
    matrix = analyzePhase('matrix', matrixMask, validMask, xCenters, yCenters, dx, dy, opt, cutCell.matrix);

    cellArea = dx * dy;
    validArea = nnz(validMask) * cellArea;
    poreArea = pore.area;
    matrixArea = matrix.area;
    globalStats = struct();
    globalStats.validArea = validArea;
    globalStats.poreArea = poreArea;
    globalStats.matrixArea = matrixArea;
    globalStats.porosity = safeDivide(poreArea, validArea);
    globalStats.matrixFraction = safeDivide(matrixArea, validArea);
    globalStats.interfaceLength = cutCell.interfaceLength;
    globalStats.specificInterfaceBulk = safeDivide(globalStats.interfaceLength, validArea);
    globalStats.specificInterfacePore = safeDivide(globalStats.interfaceLength, poreArea);
    globalStats.specificInterfaceMatrix = safeDivide(globalStats.interfaceLength, matrixArea);
    stats = pd_network_build_advanced_stats( ...
        globalStats, pore, matrix, poreMask, matrixMask, validMask, ...
        xCenters, yCenters, dx, dy, step.timestep, chunkFile, opt, cutCell);
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
