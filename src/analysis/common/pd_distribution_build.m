function data = pd_distribution_build(x, valueRange, binSize, fitTypes, emptyBinMode)
%PD_DISTRIBUTION_BUILD Build histogram, PDF, CDF, and fitted curves.
%   This is the canonical positive-valued distribution engine shared by
%   cluster and network analyses.

    if nargin < 2
        valueRange = [];
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

    x = x(:);
    x = x(isfinite(x) & (x > 0));
    data = emptyDistribution(emptyBinMode);
    if isempty(x)
        return;
    end

    edges = buildEdges(x, valueRange, binSize);
    rawCount = histcounts(x, edges);
    rawCenters = 0.5 * (edges(1:end-1) + edges(2:end));
    widths = diff(edges);
    totalCount = sum(rawCount);
    rawProbability = zeros(size(rawCount));
    rawPdf = zeros(size(rawCount));
    if totalCount > 0
        rawProbability = rawCount ./ totalCount;
        rawPdf = rawCount ./ (totalCount .* widths);
    end

    [centers, count, probability, pdfValues, keptBins] = applyEmptyBinMode( ...
        rawCenters, rawCount, rawProbability, rawPdf, emptyBinMode);
    sortedValues = sort(x, 'ascend');
    cdfCount = (1:numel(sortedValues)).';

    data.edges = edges;
    data.centers = centers;
    data.count = count;
    data.probability = probability;
    data.pdf = pdfValues;
    data.rawCenters = rawCenters;
    data.rawCount = rawCount;
    data.rawProbability = rawProbability;
    data.rawPdf = rawPdf;
    data.keptBins = keptBins;
    data.binSize = mean(widths);
    data.totalCount = totalCount;
    data.cdfX = sortedValues;
    data.cdfCount = cdfCount;
    data.cdfProbability = cdfCount ./ numel(sortedValues);
    data.fit = buildFitCurves(x, data.binSize, totalCount, fitTypes);
end

function data = emptyDistribution(emptyBinMode)
    data = struct();
    data.edges = [];
    data.centers = [];
    data.count = [];
    data.probability = [];
    data.pdf = [];
    data.rawCenters = [];
    data.rawCount = [];
    data.rawProbability = [];
    data.rawPdf = [];
    data.keptBins = [];
    data.emptyBinMode = emptyBinMode;
    data.binSize = [];
    data.totalCount = 0;
    data.cdfX = [];
    data.cdfCount = [];
    data.cdfProbability = [];
    data.fit = emptyFitStruct();
end

function [centers, count, probability, pdfValues, keptBins] = ...
        applyEmptyBinMode(rawCenters, rawCount, rawProbability, rawPdf, mode)
    keptBins = true(size(rawCount));
    centers = rawCenters;
    count = double(rawCount);
    probability = rawProbability;
    pdfValues = rawPdf;
    zeroBins = (rawCount == 0);

    switch mode
        case 'zero'
            % Keep natural zero values.
        case 'nan'
            count(zeroBins) = NaN;
            probability(zeroBins) = NaN;
            pdfValues(zeroBins) = NaN;
        case 'remove'
            keptBins = ~zeroBins;
            centers = centers(keptBins);
            count = count(keptBins);
            probability = probability(keptBins);
            pdfValues = pdfValues(keptBins);
        otherwise
            error('pd_distribution_build:InvalidEmptyBinMode', ...
                'Empty-bin mode must be zero, nan, or remove.');
    end
end

function edges = buildEdges(x, valueRange, binSize)
    if isempty(valueRange)
        xmin = min(x);
        xmax = max(x);
    else
        xmin = min(valueRange(:));
        xmax = max(valueRange(:));
    end

    if isempty(binSize)
        binSize = chooseAutoBinSize(x);
    else
        binSize = double(binSize);
    end
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

function binSize = chooseAutoBinSize(x)
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

function fit = buildFitCurves(x, binWidth, totalCount, fitTypes)
    fit = emptyFitStruct();
    if isempty(x)
        return;
    end
    xmin = min(x);
    xmax = max(x);
    if xmin == xmax
        pad = max(1e-12, abs(xmin) * 1e-6);
        fit.x = linspace(max(eps, xmin - pad), xmax + pad, 300);
    else
        fit.x = linspace(max(eps, xmin), xmax, 300);
    end

    for i = 1:numel(fitTypes)
        typeName = fitTypes{i};
        pdfValues = nan(size(fit.x));
        params = struct();
        switch typeName
            case 'powerlaw'
                [pdfValues, params] = powerlawPdf(fit.x, x);
            case 'gamma'
                [pdfValues, params] = gammaPdfMoment(fit.x, x);
            case 'lognormal'
                [pdfValues, params] = lognormalPdfMoment(fit.x, x);
        end
        fit.(typeName) = buildFitEntry( ...
            typeName, fit.x, pdfValues, params, binWidth, totalCount);
    end
end

function fit = emptyFitStruct()
    emptyEntry = struct('enabled', false, 'name', '', 'x', [], 'pdf', [], ...
        'count', [], 'probability', [], 'params', struct());
    fit = struct();
    fit.x = [];
    fit.powerlaw = emptyEntry;
    fit.gamma = emptyEntry;
    fit.lognormal = emptyEntry;
end

function entry = buildFitEntry(typeName, x, pdfValues, params, binWidth, totalCount)
    entry = struct();
    entry.enabled = any(isfinite(pdfValues) & (pdfValues >= 0));
    entry.name = fitDisplayName(typeName);
    entry.x = x;
    entry.pdf = pdfValues;
    entry.count = pdfValues .* totalCount .* binWidth;
    entry.probability = pdfValues .* binWidth;
    entry.params = params;
    entry.count(~isfinite(entry.count)) = NaN;
    entry.probability(~isfinite(entry.probability)) = NaN;
end

function [pdfValues, params] = powerlawPdf(xGrid, xSample)
    xmin = min(xSample);
    params = struct('alpha', NaN, 'xmin', xmin);
    pdfValues = nan(size(xGrid));
    if numel(xSample) < 2 || xmin <= 0
        return;
    end
    denominator = sum(log(xSample ./ xmin));
    if ~(isfinite(denominator) && denominator > 0)
        return;
    end
    alpha = 1 + numel(xSample) / denominator;
    params.alpha = alpha;
    pdfValues = zeros(size(xGrid));
    valid = xGrid >= xmin;
    pdfValues(valid) = (alpha - 1) .* (xmin .^ (alpha - 1)) .* ...
        (xGrid(valid) .^ (-alpha));
end

function [pdfValues, params] = gammaPdfMoment(xGrid, xSample)
    sampleMean = mean(xSample);
    sampleVariance = var(xSample);
    params = struct('shape', NaN, 'scale', NaN);
    pdfValues = nan(size(xGrid));
    if ~(isfinite(sampleMean) && isfinite(sampleVariance) && ...
            sampleMean > 0 && sampleVariance > 0)
        return;
    end
    shape = (sampleMean * sampleMean) / sampleVariance;
    scale = sampleVariance / sampleMean;
    params.shape = shape;
    params.scale = scale;
    pdfValues = zeros(size(xGrid));
    valid = xGrid > 0;
    pdfValues(valid) = (xGrid(valid).^(shape - 1) .* ...
        exp(-xGrid(valid) ./ scale)) ./ (gamma(shape) .* (scale .^ shape));
    pdfValues(~isfinite(pdfValues)) = NaN;
end

function [pdfValues, params] = lognormalPdfMoment(xGrid, xSample)
    logValues = log(xSample);
    mu = mean(logValues);
    sigma = std(logValues);
    params = struct('mu', mu, 'sigma', sigma);
    pdfValues = nan(size(xGrid));
    if ~(isfinite(mu) && isfinite(sigma) && sigma > 0)
        return;
    end
    pdfValues = zeros(size(xGrid));
    valid = xGrid > 0;
    pdfValues(valid) = (1 ./ (xGrid(valid) .* sigma .* sqrt(2*pi))) .* ...
        exp(-((log(xGrid(valid)) - mu) .^ 2) ./ (2 * sigma ^ 2));
    pdfValues(~isfinite(pdfValues)) = NaN;
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
