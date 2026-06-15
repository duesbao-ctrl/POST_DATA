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

    p.addParameter('Dim', 2, @isnumeric);                % 2 or 3
    p.addParameter('Dx', 1.0, @isnumeric);               % particle spacing
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
    opt.DiameterEmptyBinMode = normalizeEmptyBinMode(opt.DiameterEmptyBinMode);
    opt.DiameterPlotStyle = lower(strtrim(toChar(opt.DiameterPlotStyle)));
    opt.DiameterFitTypes = normalizeFitTypes(opt.DiameterFitTypes, 'DiameterFitTypes');
    validateMeanPower(opt.MeanPowerM, 'MeanPowerM');
    validateMeanPower(opt.MeanPowerN, 'MeanPowerN');
    validatePositiveScalarOrEmpty(opt.DiameterHistBinSize, 'DiameterHistBinSize');
    validateRangeOrEmpty(opt.DiameterRange, 'DiameterRange');
    validateEmptyBinMode(opt.DiameterEmptyBinMode);
    validateDiameterPlotStyle(opt.DiameterPlotStyle);
    validateHistScale(opt.HistScale);

    meanDefinition = struct('m', opt.MeanPowerM, 'n', opt.MeanPowerN);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex, ...
                    'ProgressMode', opt.ProgressMode};

    step = read_chunk_step_fast(clusterPath, selectorArgs{:});
    stepIdx = step.stepIndex;
    selectionInfo = buildSelectionInfo(step, opt);
    stepData = step.data;

    col = step.colIndex;
    if ~isfield(col, opt.NcountVar)
        error('cluster_postprocess:MissingNcount', 'Missing Ncount variable "%s".', opt.NcountVar);
    end

    mask = true(size(stepData, 1), 1);
    mask = applyRange(mask, stepData, col, 'c_x', opt.Range_c_x);
    mask = applyRange(mask, stepData, col, 'c_y', opt.Range_c_y);
    mask = applyRange(mask, stepData, col, 'c_z', opt.Range_c_z);
    mask = applyRange(mask, stepData, col, 'vx',  opt.Range_vx);
    mask = applyRange(mask, stepData, col, 'vy',  opt.Range_vy);
    mask = applyRange(mask, stepData, col, 'vz',  opt.Range_vz);

    dataSel = stepData(mask, :);
    ncount = dataSel(:, col.(opt.NcountVar));

    diameterAll = equivalentDiameter(ncount, opt.Dx, opt.Dim);
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
        result.totalRows = size(stepData, 1);
        result.selectedRows = 0;
        result.colIndex = col;
        result.filteredData = dataSel;
        result.diameter = [];
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

    xVarName = matlab.lang.makeValidName(opt.XVarForMean);
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
    result.totalRows = size(stepData, 1);
    result.selectedRows = size(dataSel, 1);
    result.colIndex = col;
    result.filteredData = dataSel;
    result.diameter = diameter;
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

function mask = applyRange(mask, data, col, varName, rangeVal)
    if isempty(rangeVal)
        return;
    end
    if ~isfield(col, varName)
        warning('cluster_postprocess:MissingFilterVar', ...
            'Filter variable "%s" not found. Ignore this range filter.', varName);
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

function d = equivalentDiameter(ncount, dx, dim)
    vol = ncount .* (dx ^ dim);
    if dim == 2
        d = 2 .* sqrt(vol ./ pi);
    elseif dim == 3
        d = 2 .* ((3 .* vol) ./ (4 .* pi)).^(1/3);
    else
        error('cluster_postprocess:BadDim', 'Dim must be 2 or 3.');
    end
end

function s = basicStats(x, meanPowerM, meanPowerN)
    s = struct();
    s.n = numel(x);
    s.min = min(x);
    s.max = max(x);
    s.mean = momentRatioMean(x, meanPowerM, meanPowerN);
    s.std = std(x);
    s.median = median(x);
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
    [~, fit.powerlaw] = powerlawPdf([], x);
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
            meanD(i) = momentRatioMean(d(idx), meanPowerM, meanPowerN);
        end
    end

    out = struct();
    out.edges = edges;
    out.centers = centers;
    out.meanDiameter = meanD;
    out.count = cnt;
end

function histData = buildHistogramData(d, diameterRange, binSize, fitTypes, emptyBinMode)
    if nargin < 2
        diameterRange = [];
    end
    if nargin < 3
        binSize = [];
    end
    if nargin < 4
        fitTypes = {'powerlaw', 'gamma', 'lognormal'};
    end
    if nargin < 5
        emptyBinMode = 'zero';
    end
    d = d(:);
    d = d(isfinite(d) & (d > 0));
    if isempty(d)
        histData = emptyHistogramData();
        histData.emptyBinMode = emptyBinMode;
        return;
    end

    edges = buildHistogramEdgesFromData(d, diameterRange, binSize);
    [rawCounts, edges] = histcounts(d, edges);
    rawCenters = 0.5 * (edges(1:end-1) + edges(2:end));
    widths = diff(edges);
    if isempty(widths)
        binWidth = NaN;
    else
        binWidth = mean(widths);
    end

    totalCount = sum(rawCounts);
    if totalCount > 0
        rawProb = rawCounts / totalCount;
    else
        rawProb = zeros(size(rawCounts));
    end
    [centers, counts, prob, keptBins] = applyEmptyBinMode(rawCenters, rawCounts, rawProb, emptyBinMode);

    xg = linspace(max(min(d), eps), max(d), 300);
    histData = struct();
    histData.counts = counts;
    histData.edges = edges;
    histData.centers = centers;
    histData.rawCounts = rawCounts;
    histData.rawCenters = rawCenters;
    histData.rawProbability = rawProb;
    histData.keptBins = keptBins;
    histData.emptyBinMode = emptyBinMode;
    histData.binWidth = binWidth;
    histData.totalCount = totalCount;
    histData.prob = prob;
    histData.fitX = xg;
    histData.fit = buildDiameterFitCurves(d, xg, binWidth, totalCount, fitTypes);
end

function [centers, counts, probability, keptBins] = applyEmptyBinMode(rawCenters, rawCounts, rawProbability, mode)
    keptBins = true(size(rawCounts));
    centers = rawCenters;
    counts = double(rawCounts);
    probability = rawProbability;
    zeroBins = (rawCounts == 0);

    switch mode
        case 'zero'
            % Keep the natural histogram values.
        case 'nan'
            counts(zeroBins) = NaN;
            probability(zeroBins) = NaN;
        case 'remove'
            keptBins = ~zeroBins;
            centers = centers(keptBins);
            counts = counts(keptBins);
            probability = probability(keptBins);
    end
end

function edges = buildHistogramEdgesFromData(x, diameterRange, binSize)
    x = x(:);
    x = x(isfinite(x) & (x > 0));
    if isempty(diameterRange)
        xmin = min(x);
        xmax = max(x);
    else
        xmin = min(diameterRange(:));
        xmax = max(diameterRange(:));
    end

    if ~isempty(binSize)
        binSize = double(binSize);
        if xmin == xmax
            edges = [xmin - 0.5 * binSize, xmin + 0.5 * binSize];
            return;
        end
        edges = xmin:binSize:xmax;
        if isempty(edges)
            edges = [xmin, xmin + binSize];
        elseif edges(end) < xmax
            edges(end + 1) = edges(end) + binSize;
        end
        if numel(edges) < 2
            edges = [xmin, xmin + binSize];
        end
        return;
    end

    binSize = chooseAutoDiameterBinSize(x);
    if xmin == xmax
        edges = [xmin - 0.5 * binSize, xmin + 0.5 * binSize];
        return;
    end
    edges = xmin:binSize:xmax;
    if isempty(edges)
        edges = [xmin, xmin + binSize];
    elseif edges(end) < xmax
        edges(end + 1) = edges(end) + binSize;
    end
    if numel(edges) < 2
        edges = [xmin, xmin + binSize];
    end
end

function binSize = chooseAutoDiameterBinSize(x)
    x = sort(x(:));
    n = numel(x);
    span = max(x) - min(x);
    if n < 2 || span <= 0
        binSize = max(1e-12, max(abs(x(1)), 1) * 0.1);
        return;
    end

    q25 = percentileFromSorted(x, 25);
    q75 = percentileFromSorted(x, 75);
    iqrValue = q75 - q25;
    binSize = 2 * iqrValue / (n ^ (1/3));
    if ~(isfinite(binSize) && binSize > 0)
        binSize = span / max(1, ceil(sqrt(n)));
    end
    if ~(isfinite(binSize) && binSize > 0)
        binSize = span;
    end
end

function q = percentileFromSorted(x, pct)
    n = numel(x);
    pos = 1 + (n - 1) * pct / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end

function histData = emptyHistogramData()
    histData = struct();
    histData.counts = [];
    histData.edges = [];
    histData.centers = [];
    histData.rawCounts = [];
    histData.rawCenters = [];
    histData.rawProbability = [];
    histData.keptBins = [];
    histData.emptyBinMode = 'zero';
    histData.binWidth = [];
    histData.totalCount = 0;
    histData.prob = [];
    histData.fitX = [];
    histData.fit = emptyDiameterFitStruct();
end

function fit = buildDiameterFitCurves(xSample, xGrid, binWidth, totalCount, fitTypes)
    xSample = xSample(:);
    xSample = xSample(isfinite(xSample) & (xSample > 0));
    fit = emptyDiameterFitStruct();
    fit.x = xGrid;
    if isempty(xSample) || isempty(xGrid)
        return;
    end

    for i = 1:numel(fitTypes)
        typeName = fitTypes{i};
        pdfVals = nan(size(xGrid));
        params = struct();
        switch typeName
            case 'powerlaw'
                [pdfVals, params] = powerlawPdf(xGrid, xSample);
            case 'gamma'
                [pdfVals, params] = gammaPdfMoment(xGrid, xSample);
            case 'lognormal'
                [pdfVals, params] = lognormalPdfMoment(xGrid, xSample);
        end
        fit.(typeName) = buildFitCurveEntry(typeName, xGrid, pdfVals, params, binWidth, totalCount);
    end
end

function fit = emptyDiameterFitStruct()
    emptyEntry = struct('enabled', false, 'name', '', 'x', [], 'pdf', [], ...
        'count', [], 'probability', [], 'params', struct());
    fit = struct();
    fit.x = [];
    fit.powerlaw = emptyEntry;
    fit.gamma = emptyEntry;
    fit.lognormal = emptyEntry;
end

function entry = buildFitCurveEntry(typeName, x, pdfVals, params, binWidth, totalCount)
    entry = struct();
    entry.enabled = any(isfinite(pdfVals) & (pdfVals >= 0));
    entry.name = fitDisplayName(typeName);
    entry.x = x;
    entry.pdf = pdfVals;
    entry.count = pdfVals .* totalCount .* binWidth;
    entry.probability = pdfVals .* binWidth;
    entry.params = params;
    entry.count(~isfinite(entry.count)) = NaN;
    entry.probability(~isfinite(entry.probability)) = NaN;
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

function validateRangeOrEmpty(v, name)
    if isempty(v)
        return;
    end
    if ~(isnumeric(v) && numel(v) == 2 && all(isfinite(v(:))) && v(2) >= v(1))
        error('cluster_postprocess:BadRange', ...
            '%s must be empty or a finite [min max] range with max >= min.', name);
    end
end

function mode = normalizeEmptyBinMode(value)
    mode = lower(strtrim(toChar(value)));
    switch mode
        case {'zero', '0'}
            mode = 'zero';
        case {'nan', 'na'}
            mode = 'nan';
        case {'remove', 'delete', 'drop'}
            mode = 'remove';
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

function fitTypes = normalizeFitTypes(value, name)
    if isTextScalar(value)
        textValue = lower(strtrim(toChar(value)));
        if isempty(textValue) || strcmp(textValue, 'none') || strcmp(textValue, 'off')
            fitTypes = {};
            return;
        end
        if strcmp(textValue, 'all')
            fitTypes = {'powerlaw', 'gamma', 'lognormal'};
            return;
        end
        raw = strsplit(textValue, {',', ';', '|', ' '});
    elseif iscell(value)
        raw = value;
    else
        error('cluster_postprocess:BadFitTypes', ...
            '%s must be a string or a cell array of strings.', name);
    end

    fitTypes = {};
    for i = 1:numel(raw)
        if isempty(raw{i})
            continue;
        end
        item = lower(strtrim(toChar(raw{i})));
        if isempty(item)
            continue;
        end
        if strcmp(item, 'all')
            fitTypes = {'powerlaw', 'gamma', 'lognormal'};
            return;
        elseif strcmp(item, 'none') || strcmp(item, 'off')
            fitTypes = {};
            return;
        end
        fitTypes{end + 1} = normalizeFitTypeName(item, name); %#ok<AGROW>
    end
    fitTypes = uniqueStable(fitTypes);
end

function typeName = normalizeFitTypeName(item, name)
    switch item
        case {'powerlaw', 'power-law', 'power', 'pareto', 'powerexponent', 'power-exponent'}
            typeName = 'powerlaw';
        case {'gamma', 'gam'}
            typeName = 'gamma';
        case {'lognormal', 'log-normal', 'lognorm', 'ln'}
            typeName = 'lognormal';
        otherwise
            error('cluster_postprocess:BadFitTypes', ...
                'Unsupported %s entry "%s". Use powerlaw/gamma/lognormal/all/none.', name, item);
    end
end

function out = uniqueStable(in)
    out = {};
    for i = 1:numel(in)
        if ~any(strcmp(in{i}, out))
            out{end + 1} = in{i}; %#ok<AGROW>
        end
    end
end

function [pdfVals, params] = powerlawPdf(xGrid, xSample)
    xmin = min(xSample);
    params = struct('alpha', NaN, 'xmin', xmin);
    pdfVals = nan(size(xGrid));
    if numel(xSample) < 2 || xmin <= 0
        return;
    end

    denom = sum(log(xSample ./ xmin));
    if ~(isfinite(denom) && denom > 0)
        return;
    end
    alpha = 1 + numel(xSample) / denom;
    params.alpha = alpha;
    if isempty(xGrid)
        return;
    end
    pdfVals = zeros(size(xGrid));
    valid = xGrid >= xmin;
    pdfVals(valid) = (alpha - 1) .* (xmin .^ (alpha - 1)) .* (xGrid(valid) .^ (-alpha));
end

function [pdfVals, params] = gammaPdfMoment(xGrid, xSample)
    m = mean(xSample);
    v = var(xSample);
    params = struct('shape', NaN, 'scale', NaN);
    pdfVals = nan(size(xGrid));
    if ~(isfinite(m) && isfinite(v) && m > 0 && v > 0)
        return;
    end

    shape = (m * m) / v;
    scale = v / m;
    params.shape = shape;
    params.scale = scale;
    pdfVals = zeros(size(xGrid));
    valid = xGrid > 0;
    pdfVals(valid) = (xGrid(valid).^(shape - 1) .* exp(-xGrid(valid) ./ scale)) ./ ...
        (gamma(shape) .* (scale .^ shape));
    pdfVals(~isfinite(pdfVals)) = NaN;
end

function [pdfVals, params] = lognormalPdfMoment(xGrid, xSample)
    lx = log(xSample);
    mu = mean(lx);
    sigma = std(lx);
    params = struct('mu', mu, 'sigma', sigma);
    pdfVals = nan(size(xGrid));
    if ~(isfinite(mu) && isfinite(sigma) && sigma > 0)
        return;
    end

    pdfVals = zeros(size(xGrid));
    valid = xGrid > 0;
    pdfVals(valid) = (1 ./ (xGrid(valid) .* sigma .* sqrt(2*pi))) .* ...
        exp(-((log(xGrid(valid)) - mu) .^ 2) ./ (2 * sigma ^ 2));
    pdfVals(~isfinite(pdfVals)) = NaN;
end

function name = fitDisplayName(typeName)
    switch typeName
        case 'powerlaw'
            name = 'Power-law';
        case 'gamma'
            name = 'Gamma';
        case 'lognormal'
            name = 'Lognormal';
        otherwise
            name = typeName;
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
