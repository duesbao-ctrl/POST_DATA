function out = analyze_chunk_field(chunkFile, varargin)
%ANALYZE_CHUNK_FIELD Analyze one variable on 1D/2D/3D chunk data.
%   out = ANALYZE_CHUNK_FIELD(chunkFile, ...)
%
% Core behavior:
% 1) If variable exists in raw chunk data, plot it directly.
% 2) If not, allow derived variables:
%    T, vx, vy, vz, velocity, pressure, density,
%    Sxx, Syy, Szz, Sxy, Sxz, Syz, vonMisesS,
%    gradient / grad:<var>, strainRate
%    Raw SPH fields are used directly. MD-derived density, pressure, and
%    stress use compute_temp_stress_chunk and require dV.
% 3) Otherwise throw error.

    p = inputParser;
    p.addRequired('chunkFile', @isTextScalar);
    p.addParameter('ChunkDim', 'auto', @isTextScalar); % auto, 1d, 2d, or 3d
    p.addParameter('Variable', 'c_rho', @isTextScalar);

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

    p.addParameter('dV', [], @isnumeric); % MD density/pressure/stress only
    p.addParameter('DoPlot', true, @islogical);
    p.addParameter('PlotOptions', {}, @iscell);  % passed to plot_cloud2d/plot_line1d
    p.addParameter('CoordScale', 1, @isnumeric);
    p.addParameter('CoordRangeX', [], @isnumeric);
    p.addParameter('CoordRangeY', [], @isnumeric);
    p.addParameter('CoordRangeZ', [], @isnumeric);
    p.addParameter('GradientVariable', '', @isTextScalar);
    p.addParameter('GradientSmoothLevel', 0, @isnumeric);
    p.addParameter('StrainRateVelocityComponent', 'vz', @isTextScalar);
    p.addParameter('StrainRateDensityVariable', 'c_rho', @isTextScalar);
    p.parse(chunkFile, varargin{:});
    opt = p.Results;
    chunkFile = toChar(chunkFile);
    opt.ChunkDim = toChar(opt.ChunkDim);
    opt.Variable = toChar(opt.Variable);
    opt.SelectBy = toChar(opt.SelectBy);
    opt.SlurmPath = toChar(opt.SlurmPath);
    opt.ProgressMode = toChar(opt.ProgressMode);
    opt.GradientVariable = toChar(opt.GradientVariable);
    opt.StrainRateVelocityComponent = toChar(opt.StrainRateVelocityComponent);
    opt.StrainRateDensityVariable = toChar(opt.StrainRateDensityVariable);
    validateCoordOptions(opt);
    validateGradientOptions(opt);

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
    D = step.data;
    col = step.colIndex;
    dimensionRequested = lower(strtrim(opt.ChunkDim));
    resolvedDimension = resolveChunkDimension(dimensionRequested, col);
    opt.ChunkDim = resolvedDimension;

    varReqRaw = strtrim(opt.Variable);
    varReq = matlab.lang.makeValidName(varReqRaw);

    sourceType = 'raw';
    computedPack = [];

    [x, y, coordZ, coordOk, coordMsg] = ...
        extractCoordinates(D, col, opt.ChunkDim);
    if isfield(col, varReq)
        value = D(:, col.(varReq));
        yLabel = varReqRaw;
    else
        sourceType = 'derived';
        if isGradientRequest(varReqRaw)
            ensureOneDimensionalWithCoordinates(opt, coordOk, coordMsg, varReqRaw);
            xForDerivative = x .* opt.CoordScale;
            [value, yLabel, computedPack] = computeGradientDerived(varReqRaw, D, col, ...
                chunkFile, selectorArgs, opt.dV, xForDerivative, opt);
        elseif isStrainRateRequest(varReqRaw)
            ensureOneDimensionalWithCoordinates(opt, coordOk, coordMsg, varReqRaw);
            xForDerivative = x .* opt.CoordScale;
            [value, yLabel, computedPack] = computeStrainRateDerived(D, col, ...
                chunkFile, selectorArgs, opt.dV, xForDerivative, opt);
        else
            [value, yLabel, computedPack] = computeDerived(varReqRaw, D, col, chunkFile, selectorArgs, opt.dV);
        end
    end

    if coordOk
        [x, y, coordZ, value] = ...
            applyCoordinateTransformAndRange(x, y, coordZ, value, opt);
    end

    fig = [];
    ax = [];
    plotOut = [];
    if opt.DoPlot
        if ~coordOk
            warning('analyze_chunk_field:SkipPlotMissingCoord', '%s', coordMsg);
        else
            ttl = sprintf('%s @ timestep %g', yLabel, step.timestep);
            dimTag = lower(strtrim(opt.ChunkDim));

            if strcmp(dimTag, '1d')
                validPlot = isfinite(x) & isfinite(value);
                if any(validPlot)
                    [fig, ax, plotOut] = plot_line1d(x, value, ...
                        'Title', ttl, 'XLabel', 'x', 'YLabel', yLabel, opt.PlotOptions{:});
                else
                    warning('analyze_chunk_field:SkipPlotNoFiniteData', ...
                        'No finite x/z data found. Skip 1D plotting.');
                end
            elseif strcmp(dimTag, '2d')
                validPlot = isfinite(x) & isfinite(y) & isfinite(value);
                if any(validPlot)
                    [fig, ax, plotOut] = plot_cloud2d(x, y, value, ...
                        'Title', ttl, 'XLabel', 'x', 'YLabel', 'y', 'ColorbarLabel', yLabel, opt.PlotOptions{:});
                else
                    warning('analyze_chunk_field:SkipPlotNoFiniteData', ...
                        'No finite x/y/z data found. Skip 2D plotting.');
                end
            elseif strcmp(dimTag, '3d')
                validPlot = isfinite(x) & isfinite(y) & ...
                    isfinite(coordZ) & isfinite(value);
                if any(validPlot)
                    fig = figure('Name', ttl);
                    ax = axes('Parent', fig);
                    plotOut = scatter3(ax, x(validPlot), y(validPlot), ...
                        coordZ(validPlot), 36, value(validPlot), 'filled');
                    xlabel(ax, 'x');
                    ylabel(ax, 'y');
                    zlabel(ax, 'z');
                    title(ax, ttl, 'Interpreter', 'none');
                    colorbar('peer', ax);
                    grid(ax, 'on');
                else
                    warning('analyze_chunk_field:SkipPlotNoFiniteData', ...
                        'No finite x/y/z/value data found. Skip 3D plotting.');
                end
            else
                error('analyze_chunk_field:BadChunkDim', ...
                    'ChunkDim must be ''1d'', ''2d'', or ''3d''.');
            end
        end
    end

    out = struct();
    out.filePath = chunkFile;
    out.selection = struct('SelectBy', opt.SelectBy, 'Index', opt.Index, ...
        'TimeStep', opt.TimeStep, 'Time', opt.Time);
    out.stepIndex = step.stepIndex;
    out.timestep = step.timestep;
    out.variableRequested = varReqRaw;
    out.variableUsed = yLabel;
    out.sourceType = sourceType;
    out.dimensionRequested = dimensionRequested;
    out.dimension = resolvedDimension;
    out.x = x;
    out.y = y;
    out.coordZ = coordZ;
    out.value = value;
    out.coordScale = opt.CoordScale;
    out.coordRangeX = opt.CoordRangeX;
    out.coordRangeY = opt.CoordRangeY;
    out.coordRangeZ = opt.CoordRangeZ;
    out.physicalTime = step.physicalTime;
    out.inputFormat = step.inputFormat;
    out.unitSystem = step.unitSystem;
    out.taskName = step.taskName;
    out.chunkKind = step.chunkKind;
    out.derivedSource = computedPack;
    out.figure = fig;
    out.axes = ax;
    out.plotOut = plotOut;
end

function dimension = resolveChunkDimension(requested, col)
    hasY = isfield(col, 'Coord2') || isfield(col, 'c_y');
    hasZ = isfield(col, 'Coord3') || isfield(col, 'c_z');
    switch lower(strtrim(requested))
        case 'auto'
            if hasZ
                dimension = '3d';
            elseif hasY
                dimension = '2d';
            else
                dimension = '1d';
            end
        case '1d'
            dimension = '1d';
        case '2d'
            if ~hasY
                error('analyze_chunk_field:MissingCoord2', ...
                    ['ChunkDim is 2d, but Coord2/c_y is absent. Select auto or 1d ', ...
                     'for one-dimensional chunk data.']);
            end
            dimension = '2d';
        case '3d'
            if ~hasY || ~hasZ
                error('analyze_chunk_field:MissingCoord3', ...
                    'ChunkDim is 3d, but Coord2/Coord3 coordinates are absent.');
            end
            dimension = '3d';
        otherwise
            error('analyze_chunk_field:BadChunkDim', ...
                'ChunkDim must be auto, 1d, 2d, or 3d.');
    end
end

function validateCoordOptions(opt)
    if ~(isscalar(opt.CoordScale) && isnumeric(opt.CoordScale) && isfinite(opt.CoordScale) && opt.CoordScale > 0)
        error('analyze_chunk_field:BadCoordScale', ...
            'CoordScale must be a positive finite scalar.');
    end
    validateRange(opt.CoordRangeX, 'CoordRangeX');
    validateRange(opt.CoordRangeY, 'CoordRangeY');
    validateRange(opt.CoordRangeZ, 'CoordRangeZ');
end

function validateGradientOptions(opt)
    if ~(isscalar(opt.GradientSmoothLevel) && isnumeric(opt.GradientSmoothLevel) && ...
            isfinite(opt.GradientSmoothLevel) && opt.GradientSmoothLevel >= 0)
        error('analyze_chunk_field:BadGradientSmoothLevel', ...
            'GradientSmoothLevel must be a non-negative finite scalar.');
    end
end

function validateRange(v, name)
    if isempty(v)
        return;
    end
    if ~(isnumeric(v) && numel(v) == 2 && all(isfinite(v(:))) && (v(2) >= v(1)))
        error('analyze_chunk_field:BadCoordRange', ...
            '%s must be empty or a finite [min max] range with max >= min.', name);
    end
end

function [x, y, coordZ, value] = ...
        applyCoordinateTransformAndRange(x, y, coordZ, value, opt)
    x = x .* opt.CoordScale;
    if ~isempty(y)
        y = y .* opt.CoordScale;
    end
    if ~isempty(coordZ)
        coordZ = coordZ .* opt.CoordScale;
    end

    mask = true(size(value));
    if ~isempty(opt.CoordRangeX)
        mask = mask & (x >= opt.CoordRangeX(1)) & (x <= opt.CoordRangeX(2));
    end
    if any(strcmpi(opt.ChunkDim, {'2d','3d'})) && ...
            ~isempty(opt.CoordRangeY) && ~isempty(y)
        mask = mask & (y >= opt.CoordRangeY(1)) & (y <= opt.CoordRangeY(2));
    end
    if strcmpi(opt.ChunkDim, '3d') && ...
            ~isempty(opt.CoordRangeZ) && ~isempty(coordZ)
        mask = mask & (coordZ >= opt.CoordRangeZ(1)) & ...
            (coordZ <= opt.CoordRangeZ(2));
    end

    x = x(mask);
    if ~isempty(y)
        y = y(mask);
    end
    if ~isempty(coordZ)
        coordZ = coordZ(mask);
    end
    value = value(mask);
end

function tf = isGradientRequest(varReqRaw)
    v = lower(strtrim(varReqRaw));
    tf = strcmp(v, 'gradient') || strcmp(v, 'grad') || ...
         ~isempty(regexp(v, '^(gradient|grad)\s*[:=]\s*.+$', 'once')) || ...
         ~isempty(regexp(v, '^(gradient|grad)\s*\(.+\)$', 'once'));
end

function tf = isStrainRateRequest(varReqRaw)
    v = regexprep(strtrim(varReqRaw), '[^a-zA-Z0-9]', '');
    tf = any(strcmpi(v, {'strainrate', 'strainrate1d', 'epsdot', 'epsilondot'}));
end

function ensureOneDimensionalWithCoordinates(opt, coordOk, coordMsg, varReqRaw)
    if ~strcmpi(opt.ChunkDim, '1d')
        error('analyze_chunk_field:GradientRequires1D', ...
            'Derived variable "%s" requires ChunkDim=''1d''.', varReqRaw);
    end
    if ~coordOk
        error('analyze_chunk_field:GradientMissingCoord', '%s', coordMsg);
    end
end

function targetRaw = resolveGradientTarget(varReqRaw, gradientVariable)
    targetRaw = '';
    expr = strtrim(varReqRaw);
    token = regexp(expr, '^(?:gradient|grad)\s*[:=]\s*(.+)$', 'tokens', 'once', 'ignorecase');
    if isempty(token)
        token = regexp(expr, '^(?:gradient|grad)\s*\((.+)\)$', 'tokens', 'once', 'ignorecase');
    end
    if ~isempty(token)
        targetRaw = strtrim(token{1});
    elseif ~isempty(strtrim(gradientVariable))
        targetRaw = strtrim(gradientVariable);
    end
    if isempty(targetRaw)
        error('analyze_chunk_field:MissingGradientVariable', ...
            ['Gradient calculation requires GradientVariable, or use ', ...
             'Variable=''grad:<variableName>''.']);
    end
end

function [z, yLabel, pack] = computeGradientDerived(varReqRaw, D, col, ...
    chunkFile, selectorArgs, dV, xCoord, opt)

    targetRaw = resolveGradientTarget(varReqRaw, opt.GradientVariable);
    [q, qLabel, qPack] = getChunkQuantity(targetRaw, D, col, chunkFile, selectorArgs, dV);
    [z, gradMeta] = computeSmoothedGradient1d(xCoord, q, opt.GradientSmoothLevel);

    yLabel = sprintf('d(%s)/dCoord1', qLabel);
    pack = struct();
    pack.type = 'gradient1d';
    pack.formula = 'dQ/dCoord1 at fixed timestep';
    pack.targetVariable = targetRaw;
    pack.targetLabel = qLabel;
    pack.smoothLevel = opt.GradientSmoothLevel;
    pack.coordinate = xCoord(:);
    pack.quantityRaw = q(:);
    pack.quantitySmooth = gradMeta.ySmooth;
    pack.gradient = z(:);
    pack.quantitySource = qPack;
end

function [z, yLabel, pack] = computeStrainRateDerived(D, col, ...
    chunkFile, selectorArgs, dV, xCoord, opt)

    velocityVar = strtrim(opt.StrainRateVelocityComponent);
    densityVar = strtrim(opt.StrainRateDensityVariable);
    if isempty(velocityVar)
        error('analyze_chunk_field:MissingStrainRateVelocity', ...
            'StrainRateVelocityComponent cannot be empty.');
    end
    if isempty(densityVar)
        error('analyze_chunk_field:MissingStrainRateDensity', ...
            'StrainRateDensityVariable cannot be empty.');
    end

    [u, uLabel, uPack] = getChunkQuantity(velocityVar, D, col, chunkFile, selectorArgs, dV);
    [rho, rhoLabel, rhoPack] = getChunkQuantity(densityVar, D, col, chunkFile, selectorArgs, dV);

    [du, uGradMeta] = computeSmoothedGradient1d(xCoord, u, opt.GradientSmoothLevel);
    [drho, rhoGradMeta] = computeSmoothedGradient1d(xCoord, rho, opt.GradientSmoothLevel);
    z = du + safeDivide(uGradMeta.ySmooth, rhoGradMeta.ySmooth) .* drho;

    yLabel = 'Strain Rate';
    pack = struct();
    pack.type = 'strainRate1d';
    pack.formula = 'epsilonDot = du/dCoord1 + (u/rho) * d(rho)/dCoord1';
    pack.smoothLevel = opt.GradientSmoothLevel;
    pack.coordinate = xCoord(:);
    pack.velocityVariable = velocityVar;
    pack.velocityLabel = uLabel;
    pack.velocityRaw = u(:);
    pack.velocitySmooth = uGradMeta.ySmooth;
    pack.du_dCoord1 = du(:);
    pack.velocitySource = uPack;
    pack.densityVariable = densityVar;
    pack.densityLabel = rhoLabel;
    pack.densityRaw = rho(:);
    pack.densitySmooth = rhoGradMeta.ySmooth;
    pack.drho_dCoord1 = drho(:);
    pack.densitySource = rhoPack;
    pack.strainRate = z(:);
end

function [q, qLabel, qPack] = getChunkQuantity(varReqRaw, D, col, chunkFile, selectorArgs, dV)
    varReqRaw = strtrim(varReqRaw);
    if isempty(varReqRaw)
        error('analyze_chunk_field:EmptyQuantityName', ...
            'Quantity name cannot be empty.');
    end

    qPack = struct();
    qPack.requested = varReqRaw;

    [rawIdx, rawName] = findRawColumn(col, {varReqRaw});
    if ~isempty(rawIdx)
        [q, qLabel, qPack] = buildRawQuantity(D, rawIdx, rawName, qPack, 'raw');
        return;
    end

    try
        [q, qLabel, derivedPack] = computeDerived(varReqRaw, D, col, chunkFile, selectorArgs, dV);
    catch ME
        rethrow(ME);
    end
    qPack.sourceType = 'derived';
    qPack.derivedSource = derivedPack;
end

function [q, qLabel, qPack] = buildRawQuantity(D, rawIdx, rawName, qPack, sourceType)
    q = D(:, rawIdx);
    qLabel = rawName;
    qPack.sourceType = sourceType;
    qPack.columnName = rawName;
    qPack.columnIndex = rawIdx;
end

function [idx, matchedName] = findRawColumn(col, candidates)
    idx = [];
    matchedName = '';
    fields = fieldnames(col);
    for i = 1:numel(candidates)
        name = matlab.lang.makeValidName(strtrim(candidates{i}));
        if isempty(name)
            continue;
        end
        if isfield(col, name)
            idx = col.(name);
            matchedName = name;
            return;
        end
    end

    for i = 1:numel(candidates)
        name = matlab.lang.makeValidName(strtrim(candidates{i}));
        if isempty(name)
            continue;
        end
        for j = 1:numel(fields)
            if strcmpi(fields{j}, name)
                matchedName = fields{j};
                idx = col.(matchedName);
                return;
            end
        end
    end
end

function [grad, meta] = computeSmoothedGradient1d(x, y, smoothLevel)
    x = x(:);
    y = y(:);
    if numel(x) ~= numel(y)
        error('analyze_chunk_field:GradientInputSize', ...
            'Coordinate and quantity vectors must have the same length.');
    end

    grad = nan(size(y));
    ySmooth = nan(size(y));

    coordMask = isfinite(x);
    if sum(coordMask) < 2
        error('analyze_chunk_field:TooFewGradientCoordinates', ...
            'At least two finite 1D coordinates are required for gradient calculation.');
    end

    coordIdx = find(coordMask);
    [xs, orderLocal] = sort(x(coordIdx), 'ascend');
    sortedIdx = coordIdx(orderLocal);
    if any(diff(xs) <= 0)
        error('analyze_chunk_field:DuplicateGradientCoordinate', ...
            '1D coordinates must be unique and strictly increasing after sorting.');
    end

    ys = y(sortedIdx);
    ysSmooth = gaussianSmoothNaN1d(ys, smoothLevel);
    gs = finiteDifferenceSorted1d(xs, ysSmooth);

    grad(sortedIdx) = gs;
    ySmooth(sortedIdx) = ysSmooth;

    meta = struct();
    meta.xSorted = xs;
    meta.sortedIndex = sortedIdx;
    meta.ySmoothSorted = ysSmooth;
    meta.ySmooth = ySmooth;
end

function grad = finiteDifferenceSorted1d(x, y)
    n = numel(x);
    grad = nan(size(y));

    i = 1;
    finiteY = isfinite(y);
    while i <= n
        while i <= n && ~finiteY(i)
            i = i + 1;
        end
        if i > n
            break;
        end
        j = i;
        while j <= n && finiteY(j)
            j = j + 1;
        end
        seg = i:(j-1);
        if numel(seg) == 2
            slope = (y(seg(2)) - y(seg(1))) / (x(seg(2)) - x(seg(1)));
            grad(seg) = slope;
        elseif numel(seg) > 2
            grad(seg(1)) = (y(seg(2)) - y(seg(1))) / (x(seg(2)) - x(seg(1)));
            grad(seg(end)) = (y(seg(end)) - y(seg(end-1))) / (x(seg(end)) - x(seg(end-1)));
            for k = seg(2:end-1)
                grad(k) = (y(k+1) - y(k-1)) / (x(k+1) - x(k-1));
            end
        end
        i = j;
    end
end

function ys = gaussianSmoothNaN1d(y, sigma)
    if sigma <= 0
        ys = y;
        return;
    end

    r = max(1, ceil(3 * sigma));
    t = -r:r;
    g = exp(-(t .^ 2) / (2 * sigma ^ 2));
    g = g / sum(g);

    w = double(isfinite(y));
    y0 = y;
    y0(~isfinite(y0)) = 0;

    num = conv(y0, g, 'same');
    den = conv(w, g, 'same');

    ys = num ./ den;
    ys(den <= eps) = NaN;
end

function out = safeDivide(num, den)
    out = nan(size(num));
    mask = isfinite(num) & isfinite(den) & (den ~= 0);
    out(mask) = num(mask) ./ den(mask);
end

function [z, yLabel, pack] = computeDerived(varReqRaw, D, col, chunkFile, selectorArgs, dV)
    v = lower(strtrim(varReqRaw));

    isT = strcmp(v, 't');
    isVelocity = any(strcmp(v, {'vx','vy','vz','velocity','speed'}));
    isPressure = any(strcmp(v, {'pressure','p'}));
    isDensity = any(strcmp(v, {'density','rho'}));
    isStress = any(strcmp(v, {'sxx','syy','szz','sxy','sxz','syz','vonmisess'}));

    if ~(isT || isVelocity || isPressure || isDensity || isStress)
        error('analyze_chunk_field:UnknownVariable', ...
            ['Variable "%s" not found in raw data and not in supported derived set: ', ...
             'T,vx,vy,vz,velocity,pressure,density,Sxx,Syy,Szz,Sxy,Sxz,Syz,', ...
             'vonMisesS,gradient,grad:<var>,strainRate.'], ...
             varReqRaw);
    end

    [handledNoDV, z, yLabel, pack] = computeDerivedNoDV(v, D, col);
    if handledNoDV
        return;
    end

    needsDV = isStress || isPressure || isDensity;
    if needsDV
        if isempty(dV) || ~isscalar(dV) || ~isfinite(dV) || dV <= 0
            error('analyze_chunk_field:NeedDV', ...
                'Derived variable "%s" requires a positive scalar dV.', varReqRaw);
        end
        c = compute_temp_stress_chunk(chunkFile, selectorArgs{:}, 'dV', dV);
    else
        error('analyze_chunk_field:InternalNoDVSwitch', ...
            'Unhandled dV-free derived variable: %s', varReqRaw);
    end

    switch v
        case {'pressure','p'}
            assertPressureComputed(c, 'Pressure');
            z = c.pressure;
            yLabel = 'Pressure (GPa)';

        case {'density','rho'}
            z = c.density;
            yLabel = 'Density (g/cm^3)';

        case 'sxx'
            assertComputed(c, 'xx', 'Sxx');
            z = c.sigma.xx;
            yLabel = 'Sxx';

        case 'syy'
            assertComputed(c, 'yy', 'Syy');
            z = c.sigma.yy;
            yLabel = 'Syy';

        case 'szz'
            assertComputed(c, 'zz', 'Szz');
            z = c.sigma.zz;
            yLabel = 'Szz';

        case 'sxy'
            assertComputed(c, 'xy', 'Sxy');
            z = c.sigma.xy;
            yLabel = 'Sxy';

        case 'sxz'
            assertComputed(c, 'xz', 'Sxz');
            z = c.sigma.xz;
            yLabel = 'Sxz';

        case 'syz'
            assertComputed(c, 'yz', 'Syz');
            z = c.sigma.yz;
            yLabel = 'Syz';

        case 'vonmisess'
            need = {'xx','yy','zz','xy','xz','yz'};
            for i = 1:numel(need)
                assertComputed(c, need{i}, ['vonMisesS (', need{i}, ')']);
            end
            sx = c.sigma.xx; sy = c.sigma.yy; sz = c.sigma.zz;
            txy = c.sigma.xy; txz = c.sigma.xz; tyz = c.sigma.yz;
            z = sqrt(0.5 .* ((sx-sy).^2 + (sy-sz).^2 + (sz-sx).^2 + 6.*(txy.^2 + txz.^2 + tyz.^2)));
            yLabel = 'vonMisesS';

        otherwise
            error('analyze_chunk_field:InternalVarSwitch', 'Unhandled derived variable: %s', varReqRaw);
    end

    pack = c;
end

function [handled, z, yLabel, pack] = computeDerivedNoDV(v, D, col)
    handled = true;
    pack = [];

    switch v
        case {'vx', 'vy', 'vz', 'velocity', 'speed', 't'}
            [vx, vy, vz, speed, velocityPack] = velocityFromRawColumns(D, col);
            switch v
                case 'vx'
                    z = vx;
                    yLabel = 'Vx';
                case 'vy'
                    z = vy;
                    yLabel = 'Vy';
                case 'vz'
                    z = vz;
                    yLabel = 'Vz';
                case {'velocity', 'speed'}
                    z = speed;
                    yLabel = 'Velocity';
                case 't'
                    z = temperatureFromRawColumns(D, col, vx, vy, vz);
                    yLabel = 'T';
            end
            pack = velocityPack;
            pack.sourceType = 'derived-no-dV';
            pack.requiresDV = false;

        otherwise
            handled = false;
            z = [];
            yLabel = '';
    end
end

function [vx, vy, vz, speed, pack] = velocityFromRawColumns(D, col)
    idxM = requireRawColumn(col, {'v_mass', 'mass'}, 'velocity');
    idxMvx = requireRawColumn(col, {'v_mvx', 'mvx'}, 'velocity');
    idxMvy = requireRawColumn(col, {'v_mvy', 'mvy'}, 'velocity');
    idxMvz = requireRawColumn(col, {'v_mvz', 'mvz'}, 'velocity');

    M = D(:, idxM);
    mvx = D(:, idxMvx);
    mvy = D(:, idxMvy);
    mvz = D(:, idxMvz);

    vx = safeDivide(mvx, M);
    vy = safeDivide(mvy, M);
    vz = safeDivide(mvz, M);
    speed = sqrt(vx .^ 2 + vy .^ 2 + vz .^ 2);

    pack = struct();
    pack.columns = struct('mass', idxM, 'mvx', idxMvx, 'mvy', idxMvy, 'mvz', idxMvz);
end

function T = temperatureFromRawColumns(D, col, vx, vy, vz)
    idxN = requireRawColumn(col, {'Ncount'}, 'T');
    idxM = requireRawColumn(col, {'v_mass', 'mass'}, 'T');
    idxMvx = requireRawColumn(col, {'v_mvx', 'mvx'}, 'T');
    idxMvy = requireRawColumn(col, {'v_mvy', 'mvy'}, 'T');
    idxMvz = requireRawColumn(col, {'v_mvz', 'mvz'}, 'T');
    idxKE = requireRawColumn(col, {'c_ke', 'ke'}, 'T');

    kB = 8.617343e-5;
    mvv2e = 1.0364269e-4;

    N = D(:, idxN);
    M = D(:, idxM);
    mvx = D(:, idxMvx);
    mvy = D(:, idxMvy);
    mvz = D(:, idxMvz);
    KE = D(:, idxKE);

    KE_bulk = 0.5 .* mvv2e .* M .* (vx .^ 2 + vy .^ 2 + vz .^ 2);
    KE_th = KE + KE_bulk - mvv2e .* (vx .* mvx + vy .* mvy + vz .* mvz);
    T = safeDivide(2 .* KE_th, 3 .* N .* kB);
end

function idx = requireRawColumn(col, candidates, alias)
    [idx, ~] = findRawColumn(col, candidates);
    if isempty(idx)
        error('analyze_chunk_field:MissingDerivedColumn', ...
            'Cannot compute %s without raw column(s): %s.', alias, strjoin(candidates, ', '));
    end
end

function assertComputed(c, key, alias)
    if ~isfield(c, 'computed') || ~isfield(c.computed, key) || ~c.computed.(key)
        error('analyze_chunk_field:MissingStressComponent', ...
            'Cannot compute %s because required stress component "%s" is unavailable.', alias, key);
    end
end

function assertPressureComputed(c, alias)
    if ~isfield(c, 'computed') || ~isfield(c.computed, 'pressure') || ~c.computed.pressure
        error('analyze_chunk_field:MissingPressure', ...
            'Cannot compute %s because required diagonal stress components are unavailable.', alias);
    end
end

function [x, y, coordZ, ok, msg] = extractCoordinates(D, col, chunkDim)
    ok = true;
    msg = '';
    if isfield(col, 'Coord1')
        x = D(:, col.Coord1);
    elseif isfield(col, 'c_x')
        x = D(:, col.c_x);
    elseif isfield(col, 'Chunk')
        x = D(:, col.Chunk);
    else
        x = [];
        y = [];
        coordZ = [];
        ok = false;
        msg = 'Missing x column (Coord1/c_x/Chunk). Skip plotting.';
        return;
    end

    if any(strcmpi(strtrim(chunkDim), {'2d','3d'}))
        if isfield(col, 'Coord2')
            y = D(:, col.Coord2);
        elseif isfield(col, 'c_y')
            y = D(:, col.c_y);
        else
            y = [];
            coordZ = [];
            ok = false;
            msg = '2D/3D mode requires y column (Coord2/c_y). Skip plotting.';
            return;
        end
    else
        y = [];
    end

    if strcmpi(strtrim(chunkDim), '3d')
        if isfield(col, 'Coord3')
            coordZ = D(:, col.Coord3);
        elseif isfield(col, 'c_z')
            coordZ = D(:, col.c_z);
        else
            coordZ = [];
            ok = false;
            msg = '3D mode requires z column (Coord3/c_z). Skip plotting.';
        end
    else
        coordZ = [];
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
