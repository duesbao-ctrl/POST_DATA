function stats = pd_network_build_advanced_stats(globalStats, pore, matrix, poreMask, matrixMask, validMask, ...
        xCenters, yCenters, dx, dy, timestep, chunkFile, opt, cutCell)
    if nargin < 14
        cutCell = [];
    end
    morphology = pd_network_build_morphology_stats(poreMask, matrixMask, ...
        pore.area, matrix.area, dx, dy, opt.Boundary, ...
        opt.EnableSkeletonGraph);
    toolbox = morphology.toolbox;

    stats = struct();
    stats.meta = buildPaperMeta(opt, dx, dy, timestep, chunkFile, xCenters, yCenters, toolbox, validMask, cutCell);
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
    stats.network = morphology.network;
    stats.thickness = morphology.thickness;
    stats.fragmentation = struct( ...
        'matrix', buildFragmentationStats(matrix, stats.size.matrix, stats.topology.matrix));
end

function meta = buildPaperMeta(opt, dx, dy, timestep, chunkFile, xCenters, yCenters, toolbox, validMask, cutCell)
    if nargin < 10
        cutCell = [];
    end
    meta = struct();
    meta.filePath = chunkFile;
    meta.timestep = timestep;
    meta.boundary = opt.Boundary;
    meta.thresholdN = opt.ThresholdN;
    meta.foregroundConnectivity = 4;
    meta.backgroundHoleConnectivity = 8;
    meta.meanDefinition = struct('m', opt.MeanPowerM, 'n', opt.MeanPowerN);
    meta.diameterRange = opt.DiameterRange;
    meta.diameterHistogramBinSize = opt.DiameterHistBinSize;
    meta.diameterEmptyBinMode = opt.DiameterEmptyBinMode;
    meta.diameterPlotStyle = opt.DiameterPlotStyle;
    meta.cutCell = buildCutCellMeta(opt, cutCell);
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

function meta = buildCutCellMeta(opt, cutCell)
    meta = struct();
    meta.geometryMode = opt.GeometryMode;
    meta.enabled = strcmp(opt.GeometryMode, 'cutcell');
    meta.method = opt.CutCellMethod;
    meta.fallback = opt.CutCellFallback;
    meta.plotRefinement = max(1, round(opt.CutCellPlotRefinement));
    meta.geometryFieldsAffected = ['porosity, phase area, equivalent diameter, ', ...
        'size distribution, interface length, specific interface, profiles, and phase plotting'];
    meta.topologyFieldsAffected = 'none; topology/connectivity still use the binary Ncount < ThresholdN mask.';
    meta.note = '';
    meta.fallbackCellCount = NaN;
    if isstruct(cutCell) && isfield(cutCell, 'note')
        meta.geometryMode = cutCell.geometryMode;
        meta.enabled = cutCell.enabled;
        meta.method = cutCell.method;
        meta.fallback = cutCell.fallback;
        meta.plotRefinement = cutCell.plotRefinement;
        meta.note = cutCell.note;
        meta.fallbackCellCount = cutCell.fallbackCellCount;
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

function conn = buildConnectivityPhaseStats(phase, boundary)
    conn = struct();
    conn.componentCount = phase.numComponents;
    conn.largestArea = phase.largestComponentArea;
    conn.largestFraction = phase.largestComponentFraction;
    conn.wrapsX = phase.wrapsX;
    conn.wrapsY = phase.wrapsY;
    if pd_network_is_periodic_axis(boundary, 'x')
        conn.percolatesX = logical(phase.connectivity.x.isConnected);
    else
        conn.percolatesX = logical(phase.percolatesX);
    end
    if pd_network_is_periodic_axis(boundary, 'y')
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
    [diameter, area] = pd_network_filter_component_sizes( ...
        phase.components.equivDiameter(:), phase.components.area(:), opt.DiameterRange);

    sizeStats = struct();
    sizeStats.diameterRange = opt.DiameterRange;
    sizeStats.count = numel(diameter);
    sizeStats.diameter = diameter;
    sizeStats.maxDiameter = zeroIfEmpty(maxOrNaN(diameter));
    sizeStats.meanDiameter = pd_stats_moment_ratio(diameter, opt.MeanPowerM, opt.MeanPowerN);
    sizeStats.areaWeightedMeanDiameter = weightedMean(diameter, area);
    sizeStats.hist = struct( ...
        'diameter', buildPdfHistogramData(diameter, opt.DiameterRange, ...
            opt.DiameterHistBinSize, opt.DiameterFitTypes, opt.DiameterEmptyBinMode));
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

function histData = buildPdfHistogramData(x, diameterRange, binSize, fitTypes, emptyBinMode)
    shared = pd_distribution_build(x, diameterRange, binSize, fitTypes, emptyBinMode);
    histData = struct('edges', shared.edges, 'centers', shared.centers, ...
        'count', shared.count, 'pdf', shared.pdf, ...
        'probability', shared.probability, 'rawCenters', shared.rawCenters, ...
        'rawCount', shared.rawCount, 'rawPdf', shared.rawPdf, ...
        'rawProbability', shared.rawProbability, 'keptBins', shared.keptBins, ...
        'emptyBinMode', shared.emptyBinMode, 'binSize', shared.binSize, ...
        'fit', shared.fit);
end

function frag = buildFragmentationStats(matrix, matrixSize, matrixTopology)
    frag = struct();
    frag.count = matrixTopology.beta0;
    frag.largestFraction = zeroIfNaN(matrix.largestComponentFraction);
    frag.sizeDist = matrixSize.hist.diameter;
end

function y = zeroIfNaN(x)
    if isempty(x) || (isscalar(x) && isnan(x))
        y = 0;
    else
        y = x;
    end
end
