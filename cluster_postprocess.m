function result = cluster_postprocess(clusterPath, varargin)
%CLUSTER_POSTPROCESS Cluster post-processing (selection, filtering, stats, plots).
% MATLAB R2016b compatible.
%
% Example:
% result = cluster_postprocess('cluster_chunk.txt', ...
%     'SelectBy', 'Time', 'Time', 33.0, 'SlurmPath', 'slurm_9.log', ...
%     'Dim', 2, 'Dx', 0.025, 'Range_c_x', [0 20], ...
%     'Range_diameter', [0.5 5.0], 'MeanPowerM', 1, 'MeanPowerN', 0, ...
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
    p.addParameter('Range_diameter', [], @isnumeric);

    p.addParameter('XVarForMean', 'c_x', @isTextScalar);       % c_x/c_y/c_z/vx/vy/vz
    p.addParameter('MeanNumBins', 20, @isnumeric);
    p.addParameter('MeanPowerM', 1, @isnumeric);
    p.addParameter('MeanPowerN', 0, @isnumeric);
    p.addParameter('HistNumBins', 30, @isnumeric);
    p.addParameter('DiameterHistBinSize', [], @isnumeric);
    p.addParameter('DiameterPlotRange', [], @isnumeric);
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
    validateMeanPower(opt.MeanPowerM, 'MeanPowerM');
    validateMeanPower(opt.MeanPowerN, 'MeanPowerN');
    validatePositiveScalarOrEmpty(opt.DiameterHistBinSize, 'DiameterHistBinSize');
    validateRangeOrEmpty(opt.DiameterPlotRange, 'DiameterPlotRange');

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
    diaMask = applyNumericRange(true(size(diameter)), diameter, opt.Range_diameter, 'diameter');
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
    histData = buildHistogramData(diameter, opt.HistNumBins, opt.DiameterHistBinSize);

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
        plots.countFig = plotCountDistributionWithFits(histData, fit, opt.HistScale, opt.DiameterPlotRange);
        plots.probFig = plotProbabilityDistributionWithFits(histData, fit, opt.HistScale, opt.DiameterPlotRange);
        plots.histFig = plots.countFig; % Deprecated compatibility alias.
        plots.cdfFig  = plotCDFBothDirections(diameter, opt.DiameterPlotRange);
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

function histData = buildHistogramData(d, nBins, binSize)
    if nargin < 3
        binSize = [];
    end
    edges = buildHistogramEdgesFromData(d, nBins, binSize);
    [counts, edges] = histcounts(d, edges);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    widths = diff(edges);
    if isempty(widths)
        binWidth = NaN;
    else
        binWidth = mean(widths);
    end

    totalCount = sum(counts);
    if totalCount > 0
        prob = counts / totalCount;
    else
        prob = zeros(size(counts));
    end

    xg = linspace(max(min(d), eps), max(d), 300);
    histData = struct();
    histData.counts = counts;
    histData.edges = edges;
    histData.centers = centers;
    histData.binWidth = binWidth;
    histData.totalCount = totalCount;
    histData.prob = prob;
    histData.fitX = xg;
end

function edges = buildHistogramEdgesFromData(x, nBins, binSize)
    x = x(:);
    x = x(isfinite(x) & (x > 0));
    nBins = max(1, round(nBins));
    xmin = min(x);
    xmax = max(x);

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

    if xmin == xmax
        pad = max(1e-12, abs(xmin) * 1e-6);
        edges = linspace(xmin - pad, xmax + pad, nBins + 1);
    else
        edges = linspace(xmin, xmax, nBins + 1);
    end
end

function histData = emptyHistogramData()
    histData = struct();
    histData.counts = [];
    histData.edges = [];
    histData.centers = [];
    histData.binWidth = [];
    histData.totalCount = 0;
    histData.prob = [];
    histData.fitX = [];
end

function fig = plotCountDistributionWithFits(histData, fit, histScale, xRange)
    fig = figure('Name', 'Cluster Diameter Count Distribution');

    hBar = bar(histData.centers, histData.counts, 1.0, ...
        'FaceColor', [0.35 0.6 0.85], 'EdgeColor', 'none');
    hold on;

    pLogn = lognormalPdf(histData.fitX, fit.lognormal.mu, fit.lognormal.sigma);
    pGam = gammaPdfImageForm(histData.fitX, fit.gamma.n, fit.gamma.xMean);
    scale = histData.totalCount * histData.binWidth;
    hLogn = plot(histData.fitX, pLogn * scale, 'r-', 'LineWidth', 1.6);
    hGam = plot(histData.fitX, pGam * scale, 'm--', 'LineWidth', 1.6);

    xlabel('Equivalent Diameter');
    ylabel('Count');
    title('Cluster count distribution with Lognormal/Gamma fits');
    legend([hBar, hLogn, hGam], ...
        {'Count (hist)', 'Lognormal fit', 'Gamma fit'}, ...
        'Location', 'best');
    applyHistScale(histScale);
    applyDistributionXRange(gca, xRange, histData.centers);
    grid on;
end

function fig = plotProbabilityDistributionWithFits(histData, fit, histScale, xRange)
    fig = figure('Name', 'Cluster Diameter Probability Distribution');

    hBar = bar(histData.centers, histData.prob, 1.0, ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    hold on;

    pLogn = lognormalPdf(histData.fitX, fit.lognormal.mu, fit.lognormal.sigma);
    pGam = gammaPdfImageForm(histData.fitX, fit.gamma.n, fit.gamma.xMean);
    scale = histData.binWidth;
    hLogn = plot(histData.fitX, pLogn * scale, 'r-', 'LineWidth', 1.6);
    hGam = plot(histData.fitX, pGam * scale, 'm--', 'LineWidth', 1.6);

    xlabel('Equivalent Diameter');
    ylabel('Probability');
    title('Cluster probability distribution with Lognormal/Gamma fits');
    legend([hBar, hLogn, hGam], ...
        {'Probability (bin)', 'Lognormal fit', 'Gamma fit'}, ...
        'Location', 'best');
    applyHistScale(histScale);
    applyDistributionXRange(gca, xRange, histData.centers);
    grid on;
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
    plots = struct('histFig', [], 'countFig', [], 'probFig', [], 'cdfFig', [], 'meanFig', []);
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

function p = lognormalPdf(x, mu, sigma)
    if sigma <= 0 || ~isfinite(sigma)
        p = zeros(size(x));
        return;
    end
    p = (1 ./ (x .* sigma .* sqrt(2*pi))) .* exp(-((log(x)-mu).^2) ./ (2*sigma^2));
    p(~isfinite(p)) = 0;
end

function p = gammaPdfImageForm(x, n, xMean)
%GAMMAPDFIMAGEFORM Gamma PDF in the form shown by user:
%   f(xi) = n^n/Gamma(n) * xi^(n-1) * exp(-n*xi),  xi = x/xMean
% Converted back to x-space by: f_x(x) = f_xi(x/xMean) / xMean

    if ~(isfinite(n) && isfinite(xMean)) || n <= 0 || xMean <= 0
        p = zeros(size(x));
        return;
    end

    xi = x ./ xMean;
    pXi = (n.^n ./ gamma(n)) .* (xi.^(n-1)) .* exp(-n .* xi);
    p = pXi ./ xMean;
    p(~isfinite(p)) = 0;
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
