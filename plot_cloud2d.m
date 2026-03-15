function [fig, ax, out] = plot_cloud2d(x, y, z, varargin)
%PLOT_CLOUD2D Draw a smooth 2D cloud map from x/y/z scattered or gridded data.
% MATLAB R2016b compatible.
%
% Usage:
%   [fig, ax, out] = plot_cloud2d(x, y, z, ...)
%   [fig, ax, out] = plot_cloud2d(data3col, [], [], ...)
%
% Options (Name-Value):
%   'NaNMode'        : 'color' (default) | 'transparent' | 'white'
%   'NaNColor'       : [r g b], default [0.92 0.92 0.92]
%   'SmoothLevel'    : >=0, gaussian sigma in grid cells, default 1.0
%   'PreserveNaN'    : true (default) to keep original NaN locations after smoothing
%   'CLim'           : [zmin zmax], default auto from data
%   'GridSize'       : [nx ny], default auto
%   'InterpMethod'   : 'natural' (default) | 'linear' | 'nearest'
%   'Colormap'       : colormap name/function, default 'turbo'
%   'Title'          : title text
%   'XLabel'         : x-axis label, default 'x'
%   'YLabel'         : y-axis label, default 'y'
%   'ColorbarLabel'  : colorbar label, default 'z'
%
% Output out:
%   out.X, out.Y, out.ZRaw, out.ZSmooth, out.MaskNaN

    [x, y, z] = normalizeInput(x, y, z);

    p = inputParser;
    p.addParameter('NaNMode', 'color', @isTextScalar);
    p.addParameter('NaNColor', [0.92 0.92 0.92], @(v) isnumeric(v) && numel(v)==3);
    p.addParameter('SmoothLevel', 1.0, @(v) isnumeric(v) && isscalar(v) && v>=0);
    p.addParameter('PreserveNaN', true, @islogical);
    p.addParameter('CLim', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    p.addParameter('GridSize', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    p.addParameter('InterpMethod', 'natural', @isTextScalar);
    p.addParameter('Colormap', 'turbo');
    p.addParameter('Title', '', @isTextScalar);
    p.addParameter('XLabel', 'x', @isTextScalar);
    p.addParameter('YLabel', 'y', @isTextScalar);
    p.addParameter('ColorbarLabel', 'z', @isTextScalar);
    p.parse(varargin{:});
    opt = p.Results;
    opt.NaNMode = toChar(opt.NaNMode);
    opt.InterpMethod = toChar(opt.InterpMethod);
    opt.Title = toChar(opt.Title);
    opt.XLabel = toChar(opt.XLabel);
    opt.YLabel = toChar(opt.YLabel);
    opt.ColorbarLabel = toChar(opt.ColorbarLabel);

    [Xg, Yg, Zg] = buildGrid(x, y, z, opt.GridSize, opt.InterpMethod);
    maskNaNRaw = ~isfinite(Zg);

    Zs = nanAwareGaussianSmooth(Zg, opt.SmoothLevel);
    if opt.PreserveNaN
        Zs(maskNaNRaw) = NaN;
    end
    maskNaN = ~isfinite(Zs);

    fig = figure('Color', 'w');
    ax = axes('Parent', fig);

    h = imagesc(ax, Xg(1, :), Yg(:, 1), Zs);
    set(ax, 'YDir', 'normal');
    axis(ax, 'tight');
    axis(ax, 'equal');
    box(ax, 'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.12, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);

    applyNaNMode(ax, h, maskNaN, opt.NaNMode, opt.NaNColor);

    cm = resolveColormap(opt.Colormap);
    colormap(ax, cm);

    if isempty(opt.CLim)
        v = Zs(isfinite(Zs));
        if isempty(v)
            caxis(ax, [0 1]);
        else
            vmin = min(v);
            vmax = max(v);
            if vmin == vmax
                pad = max(1e-12, abs(vmin) * 1e-6);
                caxis(ax, [vmin - pad, vmax + pad]);
            else
                caxis(ax, [vmin vmax]);
            end
        end
    else
        cmin = opt.CLim(1);
        cmax = opt.CLim(2);
        if cmin == cmax
            pad = max(1e-12, abs(cmin) * 1e-6);
            caxis(ax, [cmin - pad, cmax + pad]);
        else
            caxis(ax, [cmin cmax]);
        end
    end

    cb = colorbar(ax);
    ylabel(cb, opt.ColorbarLabel, 'FontName', 'Times New Roman', 'FontSize', 12);

    xlabel(ax, opt.XLabel, 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, opt.YLabel, 'FontName', 'Times New Roman', 'FontSize', 13);
    if ~isempty(opt.Title)
        title(ax, opt.Title, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    end

    out = struct();
    out.X = Xg;
    out.Y = Yg;
    out.ZRaw = Zg;
    out.ZSmooth = Zs;
    out.MaskNaNRaw = maskNaNRaw;
    out.MaskNaN = maskNaN;
end

function [x, y, z] = normalizeInput(x, y, z)
    if nargin >= 1 && isnumeric(x) && ~isempty(x) && size(x,2) == 3 && (isempty(y) && isempty(z))
        y = x(:,2);
        z = x(:,3);
        x = x(:,1);
    end

    if isempty(x) || isempty(y) || isempty(z)
        error('plot_cloud2d:BadInput', 'Provide x,y,z vectors, or one Nx3 matrix.');
    end

    x = x(:);
    y = y(:);
    z = z(:);

    if ~(numel(x)==numel(y) && numel(y)==numel(z))
        error('plot_cloud2d:BadInputSize', 'x,y,z must have the same length.');
    end

    % Keep NaN in z to preserve original empty-cell mask.
    % Only remove rows with invalid coordinates.
    validXY = isfinite(x) & isfinite(y);
    x = x(validXY);
    y = y(validXY);
    z = z(validXY);

    if numel(x) < 3
        error('plot_cloud2d:TooFewPoints', 'At least 3 valid points are required.');
    end
end

function [Xg, Yg, Zg] = buildGrid(x, y, z, gridSize, interpMethod)
    [x, y, z] = aggregateDuplicatePoints(x, y, z);
    ux = unique(x);
    uy = unique(y);
    isRect = (numel(ux) * numel(uy) == numel(x));

    if isempty(gridSize)
        if isRect
            nx = numel(ux);
            ny = numel(uy);
        else
            n = numel(x);
            nx = max(80, round(sqrt(n)));
            ny = nx;
        end
    else
        nx = max(10, round(gridSize(1)));
        ny = max(10, round(gridSize(2)));
    end

    if isRect && isempty(gridSize)
        [Xg, Yg] = meshgrid(ux, uy);
        % robust reshape by sorting on y then x
        A = [x, y, z];
        A = sortrows(A, [2 1]);
        Zg = reshape(A(:,3), [numel(ux), numel(uy)]).';
    else
        xv = linspace(min(x), max(x), nx);
        yv = linspace(min(y), max(y), ny);
        [Xg, Yg] = meshgrid(xv, yv);

        m = isfinite(z);
        if sum(m) < 3
            error('plot_cloud2d:TooFewFiniteZ', ...
                'At least 3 finite z values are required for scattered interpolation.');
        end
        F = scatteredInterpolant(x(m), y(m), z(m), interpMethod, 'none');
        Zg = F(Xg, Yg);
    end
end

function [xOut, yOut, zOut] = aggregateDuplicatePoints(x, y, z)
% Aggregate duplicate (x,y) points by finite-value mean.

    XY = [x, y];
    [XYu, ~, ic] = unique(XY, 'rows');
    nU = size(XYu, 1);

    if nU == numel(x)
        xOut = x;
        yOut = y;
        zOut = z;
        return;
    end

    z0 = z;
    finiteMask = isfinite(z0);
    z0(~finiteMask) = 0;

    sumZ = accumarray(ic, z0, [nU 1], @sum, 0);
    cntZ = accumarray(ic, double(finiteMask), [nU 1], @sum, 0);

    zOut = nan(nU, 1);
    hasVal = cntZ > 0;
    zOut(hasVal) = sumZ(hasVal) ./ cntZ(hasVal);

    xOut = XYu(:, 1);
    yOut = XYu(:, 2);
end

function Zs = nanAwareGaussianSmooth(Z, sigma)
    if sigma <= 0
        Zs = Z;
        return;
    end

    r = max(1, ceil(3*sigma));
    t = -r:r;
    g1 = exp(-(t.^2)/(2*sigma^2));
    g1 = g1 / sum(g1);
    K = g1' * g1;

    W = double(isfinite(Z));
    Z0 = Z;
    Z0(~isfinite(Z0)) = 0;

    num = conv2(Z0, K, 'same');
    den = conv2(W, K, 'same');

    Zs = num ./ den;
    Zs(den <= eps) = NaN;
end

function applyNaNMode(ax, h, maskNaN, nanMode, nanColor)
    mode = lower(strtrim(nanMode));
    alpha = ones(size(maskNaN));
    alpha(maskNaN) = 0;

    switch mode
        case 'transparent'
            set(ax, 'Color', 'none');
            set(h, 'AlphaData', alpha);
        case 'white'
            set(ax, 'Color', [1 1 1]);
            set(h, 'AlphaData', alpha);
        case 'color'
            set(ax, 'Color', nanColor(:).');
            set(h, 'AlphaData', alpha);
        otherwise
            error('plot_cloud2d:BadNaNMode', 'NaNMode must be color/transparent/white.');
    end
end

function cm = resolveColormap(cmapOpt)
    if isstring(cmapOpt)
        cmapOpt = char(cmapOpt);
    end
    if ischar(cmapOpt)
        switch lower(cmapOpt)
            case 'turbo'
                if exist('turbo', 'file') == 2
                    cm = turbo(256);
                else
                    cm = parula(256);
                end
            case 'parula'
                cm = parula(256);
            case 'jet'
                cm = jet(256);
            otherwise
                try
                    f = str2func(cmapOpt);
                    cm = f(256);
                catch
                    cm = parula(256);
                end
        end
    elseif isa(cmapOpt, 'function_handle')
        cm = cmapOpt(256);
    elseif isnumeric(cmapOpt) && size(cmapOpt,2) == 3
        cm = cmapOpt;
    else
        cm = parula(256);
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
