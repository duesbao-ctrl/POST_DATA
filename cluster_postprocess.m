function result = cluster_postprocess(clusterPath, varargin)
%CLUSTER_POSTPROCESS Cluster post-processing (selection, filtering, stats, plots).
% MATLAB R2016b compatible.
%
% Example:
% result = cluster_postprocess('cluster_chunk.txt', ...
%     'SelectBy', 'Time', 'Time', 33.0, 'SlurmPath', 'slurm_9.log', ...
%     'Dim', 2, 'Dx', 0.025, 'Range_c_x', [0 20], 'XVarForMean', 'c_x');

    p = inputParser;
    p.addRequired('clusterPath', @isTextScalar);
    p.addParameter('SelectBy', 'Index', @isTextScalar);        % Time | TimeStep | Index
    p.addParameter('Time', [], @isnumeric);              % physical time from slurm
    p.addParameter('TimeStep', [], @isnumeric);          % requested timestep
    p.addParameter('Index', 1, @isnumeric);              % block index in cluster file
    p.addParameter('SlurmPath', '', @isTextScalar);
    p.addParameter('SlurmModuleIndex', 1, @isnumeric);

    p.addParameter('Dim', 2, @isnumeric);                % 2 or 3
    p.addParameter('Dx', 1.0, @isnumeric);               % particle spacing
    p.addParameter('NcountVar', 'Ncount', @isTextScalar);

    p.addParameter('Range_c_x', [], @isnumeric);
    p.addParameter('Range_c_y', [], @isnumeric);
    p.addParameter('Range_c_z', [], @isnumeric);
    p.addParameter('Range_vx', [], @isnumeric);
    p.addParameter('Range_vy', [], @isnumeric);
    p.addParameter('Range_vz', [], @isnumeric);

    p.addParameter('XVarForMean', 'c_x', @isTextScalar);       % c_x/c_y/c_z/vx/vy/vz
    p.addParameter('MeanNumBins', 20, @isnumeric);
    p.addParameter('HistNumBins', 30, @isnumeric);
    p.addParameter('HistScale', 'linear', @isTextScalar);      % linear/semilogx/semilogy/loglog/log

    p.addParameter('MakePlots', true, @islogical);
    p.parse(clusterPath, varargin{:});
    opt = p.Results;
    clusterPath = toChar(clusterPath);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.NcountVar = toChar(opt.NcountVar);
    opt.XVarForMean = toChar(opt.XVarForMean);
    opt.HistScale = toChar(opt.HistScale);

    selectorArgs = {'SelectBy', opt.SelectBy, ...
                    'Index', opt.Index, ...
                    'TimeStep', opt.TimeStep, ...
                    'Time', opt.Time, ...
                    'SlurmPath', opt.SlurmPath, ...
                    'SlurmModuleIndex', opt.SlurmModuleIndex};

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
        result.meanByBin = struct('edges', [], 'centers', [], 'meanDiameter', [], 'count', []);
        result.plots = struct('histFig', [], 'cdfFig', [], 'meanFig', []);
        return;
    end

    stats = basicStats(diameter);
    fit = fitByMoments(diameter);

    xVarName = matlab.lang.makeValidName(opt.XVarForMean);
    meanBin = struct('edges', [], 'centers', [], 'meanDiameter', [], 'count', []);
    hasMeanX = isfield(col, xVarName);
    if hasMeanX
        xVals = dataSel(:, col.(xVarName));
        meanBin = binMeanDiameter(xVals, diameter, opt.MeanNumBins);
    else
        warning('cluster_postprocess:BadXVar', ...
            'XVarForMean "%s" not found. Skip mean-diameter-by-bin plot.', opt.XVarForMean);
    end

    plots = struct();
    if opt.MakePlots
        plots.histFig = plotHistogramWithFits(diameter, opt.HistNumBins, fit, opt.HistScale);
        plots.cdfFig  = plotCDFBothDirections(diameter);
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
    result.meanByBin = meanBin;
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

function s = basicStats(x)
    s = struct();
    s.n = numel(x);
    s.min = min(x);
    s.max = max(x);
    s.mean = mean(x);
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

function out = binMeanDiameter(x, d, nBins)
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
            meanD(i) = mean(d(idx));
        end
    end

    out = struct();
    out.edges = edges;
    out.centers = centers;
    out.meanDiameter = meanD;
    out.count = cnt;
end

function fig = plotHistogramWithFits(d, nBins, fit, histScale)
    fig = figure('Name', 'Cluster Diameter Histogram');

    [counts, edges] = histcounts(d, nBins);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    bw = mean(diff(edges));
    prob = counts / sum(counts);

    yyaxis left;
    hBar = bar(centers, counts, 1.0, 'FaceColor', [0.35 0.6 0.85], 'EdgeColor', 'none');
    ylabel('Count');

    yyaxis right;
    hProb = plot(centers, prob, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 4);
    hold on;

    xg = linspace(max(min(d), eps), max(d), 300);
    pLogn = lognormalPdf(xg, fit.lognormal.mu, fit.lognormal.sigma);
    pGam = gammaPdfImageForm(xg, fit.gamma.n, fit.gamma.xMean);

    hLogn = plot(xg, pLogn * bw, 'r-', 'LineWidth', 1.6);
    hGam = plot(xg, pGam * bw, 'm--', 'LineWidth', 1.6);

    ylabel('Probability');
    xlabel('Equivalent Diameter');
    title('Histogram (Count + Probability) with Lognormal/Gamma fits');
    legend([hBar, hProb, hLogn, hGam], ...
        {'Count (hist)', 'Probability (bin)', 'Lognormal fit', 'Gamma fit'}, ...
        'Location', 'best');
    applyHistScale(histScale);
    grid on;
end

function applyHistScale(histScale)
    mode = lower(strtrim(histScale));
    ax = gca;

    switch mode
        case 'linear'
            set(ax, 'XScale', 'linear');
            yyaxis left;
            set(gca, 'YScale', 'linear');
            yyaxis right;
            set(gca, 'YScale', 'linear');
        case 'semilogx'
            set(ax, 'XScale', 'log');
            yyaxis left;
            set(gca, 'YScale', 'linear');
            yyaxis right;
            set(gca, 'YScale', 'linear');
        case {'semilogy', 'semilog'}
            set(ax, 'XScale', 'linear');
            yyaxis left;
            set(gca, 'YScale', 'log');
            yyaxis right;
            set(gca, 'YScale', 'log');
        case {'loglog', 'log'}
            set(ax, 'XScale', 'log');
            yyaxis left;
            set(gca, 'YScale', 'log');
            yyaxis right;
            set(gca, 'YScale', 'log');
        otherwise
            error('cluster_postprocess:BadHistScale', ...
                'HistScale must be linear/semilogx/semilogy/loglog/log.');
    end
end

function fig = plotCDFBothDirections(d)
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
    grid on;
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
