function result = cluster_postprocess(clusterPath, varargin)
%CLUSTER_POSTPROCESS Cluster post-processing (selection, filtering, stats, plots).
% MATLAB R2016b compatible.
%
% Example:
% result = cluster_postprocess('cluster_chunk.txt', ...
%     'SelectBy', 'Time', 'Time', 33.0, 'SlurmPath', 'slurm_9.log', ...
%     'Dim', 2, 'Dx', 0.025, 'Range_c_x', [0 20], ...
%     'DiameterRange', [0.5 5.0], 'DiameterHistBinSize', 0.1, ...
%     'MeanPowerM', 1, 'MeanPowerN', 0, ...
%     'XVarForMean', 'c_x');

    p = inputParser;
    p.addRequired('clusterPath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);        % Time | TimeStep | Index
    p.addParameter('Time', [], @isnumeric);              % physical time from slurm
    p.addParameter('TimeStep', [], @isnumeric);          % requested timestep
    p.addParameter('Index', 1, @isnumeric);              % block index in cluster file
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);
    p.addParameter('ProgressMode', 'auto', @isTextScalar);
    p.addParameter('CancelCallback', @() false, @(x) isa(x, 'function_handle'));
    p.addParameter('ProgressCallback', @(fraction, message) [], ...
        @(x) isa(x, 'function_handle'));

    p.addParameter('Dim', 2, @isnumeric);                % 2 or 3
    p.addParameter('Dx', 1.0, @isnumeric);               % particle spacing
    p.addParameter('ParticleVolume', [], @isnumeric);
    p.addParameter('ThinDirectionThickness', [], @isnumeric);
    p.addParameter('NcountVar', 'Ncount', @isTextScalar);

    p.addParameter('Range_c_x', [], @isnumeric);
    p.addParameter('Range_c_y', [], @isnumeric);
    p.addParameter('Range_c_z', [], @isnumeric);
    p.addParameter('Range_vx', [], @isnumeric);
    p.addParameter('Range_vy', [], @isnumeric);
    p.addParameter('Range_vz', [], @isnumeric);
    p.addParameter('DiameterRange', [], @isnumeric);

    p.addParameter('XVarForMean', 'c_x', @isTextScalar);       % c_x/c_y/c_z/vx/vy/vz
    p.addParameter('MeanNumBins', 20, @isnumeric);
    p.addParameter('MeanPowerM', 1, @isnumeric);
    p.addParameter('MeanPowerN', 0, @isnumeric);
    p.addParameter('DiameterHistBinSize', [], @isnumeric);
    p.addParameter('DiameterEmptyBinMode', 'zero', @isTextScalar);
    p.addParameter('DiameterPlotStyle', 'bar', @isTextScalar);
    p.addParameter('DiameterFitTypes', {'powerlaw', 'gamma', 'lognormal'});
    p.addParameter('HistScale', 'linear', @isTextScalar);      % linear/semilogx/semilogy/loglog/log

    p.addParameter('MakePlots', true, @islogical);
    p.parse(clusterPath, varargin{:});
    opt = p.Results;
    clusterPath = toChar(clusterPath);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);
    opt.NcountVar = toChar(opt.NcountVar);
    opt.XVarForMean = toChar(opt.XVarForMean);
    opt.HistScale = toChar(opt.HistScale);
    opt.DiameterEmptyBinMode = pd_distribution_normalize_empty_bin_mode( ...
        opt.DiameterEmptyBinMode);
    opt.DiameterPlotStyle = lower(strtrim(toChar(opt.DiameterPlotStyle)));
    opt.DiameterFitTypes = pd_distribution_normalize_fit_types( ...
        opt.DiameterFitTypes, 'DiameterFitTypes', ...
        'cluster_postprocess:BadFitTypes');
    validateMeanPower(opt.MeanPowerM, 'MeanPowerM');
    validateMeanPower(opt.MeanPowerN, 'MeanPowerN');
    validatePositiveScalarOrEmpty(opt.DiameterHistBinSize, 'DiameterHistBinSize');
    validateRangeOrEmpty(opt.DiameterRange, 'DiameterRange');
    validateEmptyBinMode(opt.DiameterEmptyBinMode);
    validateDiameterPlotStyle(opt.DiameterPlotStyle);
    validateHistScale(opt.HistScale);

    sizeModel = buildClusterSizeModel(opt);
    meanDefinition = struct('m', opt.MeanPowerM, 'n', opt.MeanPowerN);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex, ...
                    'ProgressMode', opt.ProgressMode, ...
                    'CancelCallback', opt.CancelCallback, ...
                    'ProgressCallback', opt.ProgressCallback};

    step = read_chunk_step_fast(clusterPath, selectorArgs{:});
    stepIdx = step.stepIndex;
    selectionInfo = buildSelectionInfo(step, opt);
    stepData = step.data;

    col = step.colIndex;
    if ~isfield(col, opt.NcountVar)
        error('cluster_postprocess:MissingNcount', 'Missing Ncount variable "%s".', opt.NcountVar);
    end

    mask = true(size(stepData, 1), 1);
    mask = applyRange(mask, stepData, col, {'c_x','x'}, opt.Range_c_x);
    mask = applyRange(mask, stepData, col, {'c_y','y'}, opt.Range_c_y);
    mask = applyRange(mask, stepData, col, {'c_z','z'}, opt.Range_c_z);
    mask = applyRange(mask, stepData, col, {'vx'}, opt.Range_vx);
    mask = applyRange(mask, stepData, col, {'vy'}, opt.Range_vy);
    mask = applyRange(mask, stepData, col, {'vz'}, opt.Range_vz);

    dataSel = stepData(mask, :);
    ncount = dataSel(:, col.(opt.NcountVar));

    diameterAll = equivalentDiameter(ncount, sizeModel);
    validDia = isfinite(diameterAll) & (diameterAll > 0);
    dataSel = dataSel(validDia, :);
    diameter = diameterAll(validDia);
    diaMask = applyNumericRange(true(size(diameter)), diameter, opt.DiameterRange, 'diameter');
    dataSel = dataSel(diaMask, :);
    diameter = diameter(diaMask);

    if isempty(diameter)
        warning('cluster_postprocess:NoValidClusters', ...
            'No valid clusters after filtering. Skip plotting.');
        result = struct();
        result.clusterPath = clusterPath;
        result.selection = selectionInfo;
        result.timestep = step.timestep;
        result.stepIndex = stepIdx;
        result.physicalTime = step.physicalTime;
        result.inputFormat = step.inputFormat;
        result.unitSystem = step.unitSystem;
        result.taskName = step.taskName;
        result.chunkKind = step.chunkKind;
        result.totalRows = size(stepData, 1);
        result.selectedRows = 0;
        result.colIndex = col;
        result.filteredData = dataSel;
        result.diameter = [];
        result.sizeModel = sizeModel;
        result.diameterRange = opt.DiameterRange;
        result.diameterEmptyBinMode = opt.DiameterEmptyBinMode;
        result.diameterPlotStyle = opt.DiameterPlotStyle;
        result.stats = [];
        result.fit = [];
        result.meanDefinition = meanDefinition;
        result.meanByBin = struct('edges', [], 'centers', [], 'meanDiameter', [], 'count', []);
        result.hist = emptyHistogramData();
        result.plots = emptyPlotsStruct();
        return;
    end

    stats = basicStats(diameter, opt.MeanPowerM, opt.MeanPowerN);
    fit = fitByMoments(diameter);
    histData = buildHistogramData(diameter, opt.DiameterRange, ...
        opt.DiameterHistBinSize, opt.DiameterFitTypes, opt.DiameterEmptyBinMode);

    xVarName = resolveClusterVariable(col, opt.XVarForMean);
    meanBin = struct('edges', [], 'centers', [], 'meanDiameter', [], 'count', []);
    hasMeanX = isfield(col, xVarName);
    if hasMeanX
        xVals = dataSel(:, col.(xVarName));
        meanBin = binMeanDiameter(xVals, diameter, opt.MeanNumBins, opt.MeanPowerM, opt.MeanPowerN);
    else
        warning('cluster_postprocess:BadXVar', ...
            'XVarForMean "%s" not found. Skip mean-diameter-by-bin plot.', opt.XVarForMean);
    end

    plots = emptyPlotsStruct();
    if opt.MakePlots
        plots.countFig = plotCountDistributionWithFits(histData, opt.HistScale, ...
            opt.DiameterRange, opt.DiameterPlotStyle);
        plots.probFig = plotProbabilityDistributionWithFits(histData, opt.HistScale, ...
            opt.DiameterRange, opt.DiameterPlotStyle);
        plots.cdfFig  = plotCDFBothDirections(diameter, opt.DiameterRange);
        hasMeanData = any(~isnan(meanBin.meanDiameter) & (meanBin.count > 0));
        if hasMeanX && hasMeanData
            plots.meanFig = plotMeanByBin(meanBin, opt.XVarForMean);
        else
            if hasMeanX && ~hasMeanData
                warning('cluster_postprocess:NoMeanPlotData', ...
                    'No valid binned mean data for "%s". Skip mean-diameter-by-bin plot.', opt.XVarForMean);
            end
            plots.meanFig = [];
        end
    end

    result = struct();
    result.clusterPath = clusterPath;
    result.selection = selectionInfo;
    result.timestep = step.timestep;
    result.stepIndex = stepIdx;
    result.physicalTime = step.physicalTime;
    result.inputFormat = step.inputFormat;
    result.unitSystem = step.unitSystem;
    result.taskName = step.taskName;
    result.chunkKind = step.chunkKind;
    result.totalRows = size(stepData, 1);
    result.selectedRows = size(dataSel, 1);
    result.colIndex = col;
    result.filteredData = dataSel;
    result.diameter = diameter;
    result.sizeModel = sizeModel;
    result.diameterRange = opt.DiameterRange;
    result.diameterEmptyBinMode = opt.DiameterEmptyBinMode;
    result.diameterPlotStyle = opt.DiameterPlotStyle;
    result.stats = stats;
    result.fit = fit;
    result.meanDefinition = meanDefinition;
    result.meanByBin = meanBin;
    result.hist = histData;
    result.plots = plots;
end

function info = buildSelectionInfo(step, opt)
    mode = lower(strtrim(opt.SelectBy));
    info = struct();
    info.mode = mode;
    info.stepIndex = step.stepIndex;
    info.mappedClusterTimeStep = step.timestep;
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

function mask = applyRange(mask, data, col, candidates, rangeVal)
    if isempty(rangeVal)
        return;
    end
    varName = resolveFirstColumn(col, candidates);
    if isempty(varName)
        warning('cluster_postprocess:MissingFilterVar', ...
            'Filter variable "%s" not found. Ignore this range filter.', ...
            strjoin(candidates, '/'));
        return;
    end
    if numel(rangeVal) ~= 2
        error('cluster_postprocess:BadRange', 'Range for %s must be [min max].', varName);
    end
    lo = min(rangeVal(:));
    hi = max(rangeVal(:));
    v = data(:, col.(varName));
    mask = mask & (v >= lo) & (v <= hi);
end

function mask = applyNumericRange(mask, values, rangeVal, varName)
    if isempty(rangeVal)
        return;
    end
    if numel(rangeVal) ~= 2
        error('cluster_postprocess:BadRange', 'Range for %s must be [min max].', varName);
    end
    lo = min(rangeVal(:));
    hi = max(rangeVal(:));
    mask = mask & (values >= lo) & (values <= hi);
end

function model = buildClusterSizeModel(opt)
    dim = opt.Dim;
    if ~(isscalar(dim) && isfinite(dim) && any(dim == [2, 3]))
        error('cluster_postprocess:BadDim', 'Dim must be 2 or 3.');
    end
    dx = opt.Dx;
    validatePositivePhysicalScalarOrEmpty(dx, 'Dx');
    particleVolume = opt.ParticleVolume;
    thinThickness = opt.ThinDirectionThickness;
    validatePositivePhysicalScalarOrEmpty(particleVolume, 'ParticleVolume');
    validatePositivePhysicalScalarOrEmpty(thinThickness, 'ThinDirectionThickness');
    if ~isempty(particleVolume) && dim ~= 3
        error('cluster_postprocess:ParticleVolumeRequires3D', ...
            'ParticleVolume is only defined for Dim=3 cluster sizing.');
    end
    if ~isempty(thinThickness) && dim ~= 3
        error('cluster_postprocess:ThinProjectionRequires3D', ...
            'ThinDirectionThickness requires Dim=3.');
    end
    if ~isempty(thinThickness) && isempty(particleVolume)
        error('cluster_postprocess:ThinProjectionRequiresParticleVolume', ...
            'ThinDirectionThickness requires ParticleVolume.');
    end
    if isempty(particleVolume) && isempty(dx)
        error('cluster_postprocess:DxRequiredForSimpleCubicSizing', ...
            'Dx is required when ParticleVolume is not specified.');
    end
    model = struct('dim', dim, 'dx', valueOrNaN(dx), 'particleVolume', ...
        valueOrNaN(particleVolume), 'thinDirectionThickness', ...
        valueOrNaN(thinThickness), 'volumePerParticle', NaN, ...
        'areaPerParticle', NaN, 'projectedAreaPerParticle', NaN, ...
        'usesParticleVolume', ~isempty(particleVolume), ...
        'usesThinProjection', ~isempty(thinThickness));
    if dim == 2
        model.mode = 'simple-cubic-2d';
        model.areaPerParticle = dx ^ 2;
        model.projectedAreaPerParticle = model.areaPerParticle;
        model.diameterDefinition = 'area=Ncount*Dx^2';
    elseif isempty(particleVolume)
        model.mode = 'simple-cubic-3d';
        model.volumePerParticle = dx ^ 3;
        model.diameterDefinition = 'volume=Ncount*Dx^3';
    elseif isempty(thinThickness)
        model.mode = 'explicit-volume-3d';
        model.volumePerParticle = particleVolume;
        model.diameterDefinition = 'volume=Ncount*ParticleVolume';
    else
        model.mode = 'quasi-2d-projection';
        model.volumePerParticle = particleVolume;
        model.projectedAreaPerParticle = particleVolume / thinThickness;
        model.diameterDefinition = 'area=Ncount*ParticleVolume/ThinDirectionThickness';
    end
end

function name = resolveClusterVariable(col, requested)
    requested = matlab.lang.makeValidName(toChar(requested));
    candidates = {requested};
    switch requested
        case 'c_x'
            candidates = {'c_x','x'};
        case 'c_y'
            candidates = {'c_y','y'};
        case 'c_z'
            candidates = {'c_z','z'};
    end
    name = resolveFirstColumn(col, candidates);
end

function name = resolveFirstColumn(col, candidates)
    name = '';
    for i = 1:numel(candidates)
        if isfield(col, candidates{i})
            name = candidates{i};
            return;
        end
    end
end

function d = equivalentDiameter(ncount, model)
    if any(strcmp(model.mode, {'simple-cubic-2d','quasi-2d-projection'}))
        if strcmp(model.mode, 'simple-cubic-2d')
            area = ncount .* model.areaPerParticle;
        else
            area = ncount .* model.projectedAreaPerParticle;
        end
        d = 2 .* sqrt(area ./ pi);
    else
        volume = ncount .* model.volumePerParticle;
        d = 2 .* ((3 .* volume) ./ (4 .* pi)).^(1/3);
    end
end

function value = valueOrNaN(value)
    if isempty(value), value = NaN; end
end

function s = basicStats(x, meanPowerM, meanPowerN)
    shared = pd_stats_summary(x, meanPowerM, meanPowerN);
    s = rmfield(shared, {'quantileProb','quantileValue'});
end

function fit = fitByMoments(x)
    lx = log(x);
    mu = mean(lx);
    sigma = std(lx);

    m = mean(x);
    v = var(x);
    if v <= 0
        n = NaN;
    else
        n = (m * m) / v;
    end

    fit = struct();
    fit.lognormal = struct('mu', mu, 'sigma', sigma);
    fit.gamma = struct('n', n, 'xMean', m);
    distribution = pd_distribution_build(x, [], [], {'powerlaw'}, 'zero');
    fit.powerlaw = distribution.fit.powerlaw.params;
end

function out = binMeanDiameter(x, d, nBins, meanPowerM, meanPowerN)
    if isempty(x)
        out = struct('edges', [], 'centers', [], 'meanDiameter', [], 'count', []);
        return;
    end

    nBins = max(1, round(nBins));
    xmin = min(x);
    xmax = max(x);
    if xmin == xmax
        edges = [xmin-0.5, xmax+0.5];
    else
        edges = linspace(xmin, xmax, nBins+1);
    end

    binId = discretize(x, edges);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));

    meanD = nan(1, numel(centers));
    cnt = zeros(1, numel(centers));
    for i = 1:numel(centers)
        idx = (binId == i);
        cnt(i) = sum(idx);
        if cnt(i) > 0
            meanD(i) = pd_stats_moment_ratio(d(idx), meanPowerM, meanPowerN);
        end
    end

    out = struct();
    out.edges = edges;
    out.centers = centers;
    out.meanDiameter = meanD;
    out.count = cnt;
end

function histData = buildHistogramData(d, diameterRange, binSize, fitTypes, emptyBinMode)
    shared = pd_distribution_build(d, diameterRange, binSize, fitTypes, emptyBinMode);
    histData = struct();
    histData.counts = shared.count;
    histData.edges = shared.edges;
    histData.centers = shared.centers;
    histData.rawCounts = shared.rawCount;
    histData.rawCenters = shared.rawCenters;
    histData.rawProbability = shared.rawProbability;
    histData.keptBins = shared.keptBins;
    histData.emptyBinMode = shared.emptyBinMode;
    histData.binWidth = shared.binSize;
    histData.totalCount = shared.totalCount;
    histData.prob = shared.probability;
    histData.fitX = shared.fit.x;
    histData.fit = shared.fit;
end

function histData = emptyHistogramData()
    histData = buildHistogramData([], [], [], {}, 'zero');
end

function [handles, labels] = plotFitCurves(ax, fit, yField)
    handles = [];
    labels = {};
    if isempty(fit) || isempty(fit.x)
        return;
    end

    types = {'powerlaw', 'gamma', 'lognormal'};
    styles = {'k-', 'm--', 'r-'};
    for i = 1:numel(types)
        entry = fit.(types{i});
        if ~entry.enabled || ~isfield(entry, yField)
            continue;
        end
        y = entry.(yField);
        valid = isfinite(entry.x) & isfinite(y) & (y >= 0);
        if ~any(valid)
            continue;
        end
        handles(end + 1) = plot(ax, entry.x(valid), y(valid), styles{i}, 'LineWidth', 1.6); %#ok<AGROW>
        labels{end + 1} = [entry.name, ' fit']; %#ok<AGROW>
    end
end

function fig = plotCountDistributionWithFits(histData, histScale, xRange, plotStyle)
    fig = figure('Name', 'Cluster Diameter Count Distribution');

    h = [];
    labels = {};
    ax = gca;
    hData = plotDistributionSeries(ax, histData.centers, histData.counts, ...
        plotStyle, [0.35 0.6 0.85]);
    hold on;
    h(end + 1) = hData;
    labels{end + 1} = 'Count (hist)';
    [hFit, fitLabels] = plotFitCurves(ax, histData.fit, 'count');
    h = [h, hFit];
    labels = [labels, fitLabels];

    xlabel('Equivalent Diameter');
    ylabel('Count');
    title('Cluster count distribution with distribution fits');
    if numel(h) > 1
        legend(h, labels, 'Location', 'best');
    end
    applyHistScale(histScale);
    applyDistributionXRange(ax, xRange, histData.centers);
    grid on;
end

function fig = plotProbabilityDistributionWithFits(histData, histScale, xRange, plotStyle)
    fig = figure('Name', 'Cluster Diameter Probability Distribution');

    h = [];
    labels = {};
    ax = gca;
    hData = plotDistributionSeries(ax, histData.centers, histData.prob, ...
        plotStyle, [0.7 0.7 0.7]);
    hold on;
    h(end + 1) = hData;
    labels{end + 1} = 'Probability (bin)';
    [hFit, fitLabels] = plotFitCurves(ax, histData.fit, 'probability');
    h = [h, hFit];
    labels = [labels, fitLabels];

    xlabel('Equivalent Diameter');
    ylabel('Probability');
    title('Cluster probability distribution with distribution fits');
    if numel(h) > 1
        legend(h, labels, 'Location', 'best');
    end
    applyHistScale(histScale);
    applyDistributionXRange(ax, xRange, histData.centers);
    grid on;
end

function h = plotDistributionSeries(ax, x, y, plotStyle, color)
    x = x(:);
    y = y(:);
    if isempty(x) || isempty(y)
        h = [];
        return;
    end

    switch plotStyle
        case 'bar'
            h = bar(ax, x, y, 1.0, 'FaceColor', color, 'EdgeColor', 'none');
        case 'scatter'
            valid = isfinite(x) & isfinite(y);
            h = scatter(ax, x(valid), y(valid), 36, color, 'filled');
    end
end

function applyHistScale(histScale)
    mode = lower(strtrim(histScale));
    ax = gca;

    switch mode
        case 'linear'
            set(ax, 'XScale', 'linear');
            set(ax, 'YScale', 'linear');
        case 'semilogx'
            set(ax, 'XScale', 'log');
            set(ax, 'YScale', 'linear');
        case {'semilogy', 'semilog'}
            set(ax, 'XScale', 'linear');
            set(ax, 'YScale', 'log');
        case {'loglog', 'log'}
            set(ax, 'XScale', 'log');
            set(ax, 'YScale', 'log');
        otherwise
            error('cluster_postprocess:BadHistScale', ...
                'HistScale must be linear/semilogx/semilogy/loglog/log.');
    end
end

function fig = plotCDFBothDirections(d, xRange)
    fig = figure('Name', 'Cluster Diameter CDF');

    xAsc = sort(d(:), 'ascend');
    n = numel(xAsc);
    countAsc = (1:n)';
    probAsc = countAsc / n;

    xDesc = flipud(xAsc);
    countDesc = (1:n)';
    probDesc = countDesc / n;

    yyaxis left;
    plot(xAsc, countAsc, 'b-', 'LineWidth', 1.4); hold on;
    plot(xDesc, countDesc, 'b--', 'LineWidth', 1.4);
    ylabel('Cumulative Count');

    yyaxis right;
    plot(xAsc, probAsc, 'r-', 'LineWidth', 1.4);
    plot(xDesc, probDesc, 'r--', 'LineWidth', 1.4);
    ylabel('Cumulative Probability');

    xlabel('Equivalent Diameter');
    title('CDF (small to large / large to small)');
    legend({'Count asc', 'Count desc', 'Prob asc', 'Prob desc'}, 'Location', 'best');
    applyDistributionXRange(gca, xRange, xAsc);
    grid on;
end

function applyDistributionXRange(ax, xRange, xFallback)
    if isempty(xRange)
        return;
    end
    xlim(ax, normalizePlotRange(xRange, xFallback));
end

function rangeOut = normalizePlotRange(rangeIn, values)
    rangeOut = double(rangeIn(:).');
    if isempty(rangeOut) || (rangeOut(2) > rangeOut(1))
        return;
    end
    values = values(:);
    values = values(isfinite(values));
    if numel(values) >= 2
        d = diff(sort(values));
        d = d(d > 0);
        if ~isempty(d)
            pad = 0.5 * min(d);
        else
            pad = 0.5;
        end
    else
        pad = 0.5;
    end
    rangeOut = [rangeOut(1) - pad, rangeOut(2) + pad];
end

function fig = plotMeanByBin(meanBin, xVarName)
    fig = figure('Name', 'Mean Diameter by Bins');

    valid = ~isnan(meanBin.meanDiameter) & (meanBin.count > 0);
    plot(meanBin.centers(valid), meanBin.meanDiameter(valid), 'o-', 'LineWidth', 1.5);
    xlabel(xVarName);
    ylabel('Mean Equivalent Diameter');
    title('Mean diameter vs binned variable');
    grid on;
end

function plots = emptyPlotsStruct()
    plots = struct('countFig', [], 'probFig', [], 'cdfFig', [], 'meanFig', []);
end

function validateMeanPower(v, name)
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v))
        error('cluster_postprocess:BadMeanPower', ...
            '%s must be a finite real scalar.', name);
    end
end

function validatePositiveScalarOrEmpty(v, name)
    if isempty(v)
        return;
    end
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0)
        error('cluster_postprocess:BadHistogramBinSize', ...
            '%s must be empty or a positive finite scalar.', name);
    end
end

function validatePositivePhysicalScalarOrEmpty(v, name)
    if isempty(v), return; end
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0)
        error('cluster_postprocess:BadPhysicalScale', ...
            '%s must be empty or a positive finite scalar.', name);
    end
end

function validateRangeOrEmpty(v, name)
    if isempty(v)
        return;
    end
    if ~(isnumeric(v) && numel(v) == 2 && all(isfinite(v(:))) && v(2) >= v(1))
        error('cluster_postprocess:BadRange', ...
            '%s must be empty or a finite [min max] range with max >= min.', name);
    end
end

function validateEmptyBinMode(mode)
    valid = {'zero', 'nan', 'remove'};
    if ~any(strcmp(mode, valid))
        error('cluster_postprocess:BadEmptyBinMode', ...
            'DiameterEmptyBinMode must be zero/nan/remove.');
    end
end

function validateDiameterPlotStyle(plotStyle)
    valid = {'bar', 'scatter'};
    if ~any(strcmp(plotStyle, valid))
        error('cluster_postprocess:BadDiameterPlotStyle', ...
            'DiameterPlotStyle must be bar or scatter.');
    end
end

function validateHistScale(histScale)
    valid = {'linear', 'semilogx', 'semilogy', 'loglog', 'semilog', 'log'};
    if ~any(strcmpi(strtrim(histScale), valid))
        error('cluster_postprocess:BadHistScale', ...
            'HistScale must be linear/semilogx/semilogy/loglog.');
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
