function [fig, ax, out] = plot_line1d(x, y, varargin)
%PLOT_LINE1D Draw a styled 1D line plot with NaN handling and smoothing.
% MATLAB R2016b compatible.
%
% Usage:
%   [fig, ax, out] = plot_line1d(x, y, ...)
%   [fig, ax, out] = plot_line1d(data2col, [], ...)
%
% Options (Name-Value):
%   'NaNMode'       : 'break' (default) | 'omit' | 'interp'
%   'SmoothLevel'   : >=0, gaussian sigma in sample points, default 0
%   'PreserveNaN'   : true (default), keep NaN locations after smoothing
%   'YLim'          : [ymin ymax], default auto
%   'XLim'          : [xmin xmax], default auto
%   'LineWidth'     : default 1.8
%   'LineColor'     : default [0.10 0.35 0.75]
%   'ShowMarkers'   : true/false, default false
%   'MarkerSize'    : default 5
%   'Title'         : title text
%   'XLabel'        : x-axis label, default 'x'
%   'YLabel'        : y-axis label, default 'y'
%
% Output out:
%   out.x, out.yRaw, out.yPlot, out.validMask

    [x, y] = normalizeInput(x, y);

    p = inputParser;
    p.addParameter('NaNMode', 'break', @isTextScalar);
    p.addParameter('SmoothLevel', 0, @(v) isnumeric(v) && isscalar(v) && v>=0);
    p.addParameter('PreserveNaN', true, @islogical);
    p.addParameter('YLim', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    p.addParameter('XLim', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    p.addParameter('LineWidth', 1.8, @(v) isnumeric(v) && isscalar(v) && v>0);
    p.addParameter('LineColor', [0.10 0.35 0.75], @(v) isnumeric(v) && numel(v)==3);
    p.addParameter('ShowMarkers', false, @islogical);
    p.addParameter('MarkerSize', 5, @(v) isnumeric(v) && isscalar(v) && v>0);
    p.addParameter('Title', '', @isTextScalar);
    p.addParameter('XLabel', 'x', @isTextScalar);
    p.addParameter('YLabel', 'y', @isTextScalar);
    p.parse(varargin{:});
    opt = p.Results;
    opt.NaNMode = toChar(opt.NaNMode);
    opt.Title = toChar(opt.Title);
    opt.XLabel = toChar(opt.XLabel);
    opt.YLabel = toChar(opt.YLabel);

    [x, idx] = sort(x, 'ascend');
    y = y(idx);
    yRaw = y;
    maskNaNRaw = ~isfinite(yRaw);

    % Always smooth with raw NaN excluded from kernel contribution.
    ySmooth = gaussianSmoothNaN1d(yRaw, opt.SmoothLevel);
    if opt.PreserveNaN
        ySmooth(maskNaNRaw) = NaN;
    end
    y = applyNaNMode(ySmooth, x, maskNaNRaw, opt.NaNMode);

    fig = figure('Color', 'w');
    ax = axes('Parent', fig);

    if opt.ShowMarkers
        plot(ax, x, y, '-o', ...
            'LineWidth', opt.LineWidth, ...
            'Color', opt.LineColor(:).', ...
            'MarkerSize', opt.MarkerSize, ...
            'MarkerFaceColor', opt.LineColor(:).');
    else
        plot(ax, x, y, '-', ...
            'LineWidth', opt.LineWidth, ...
            'Color', opt.LineColor(:).');
    end

    box(ax, 'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.16, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);

    xlabel(ax, opt.XLabel, 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, opt.YLabel, 'FontName', 'Times New Roman', 'FontSize', 13);
    if ~isempty(opt.Title)
        title(ax, opt.Title, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    end

    if ~isempty(opt.XLim)
        xlim(ax, [opt.XLim(1), opt.XLim(2)]);
    end
    if ~isempty(opt.YLim)
        ylim(ax, [opt.YLim(1), opt.YLim(2)]);
    end

    out = struct();
    out.x = x;
    out.yRaw = yRaw;
    out.yPlot = y;
    out.maskNaNRaw = maskNaNRaw;
    out.validMask = isfinite(yRaw);
end

function [x, y] = normalizeInput(x, y)
    if nargin >= 1 && isnumeric(x) && ~isempty(x) && size(x,2) == 2 && isempty(y)
        y = x(:,2);
        x = x(:,1);
    end

    if isempty(x) || isempty(y)
        error('plot_line1d:BadInput', 'Provide x,y vectors, or one Nx2 matrix.');
    end

    x = x(:);
    y = y(:);

    if numel(x) ~= numel(y)
        error('plot_line1d:BadInputSize', 'x and y must have the same length.');
    end

    goodX = isfinite(x);
    x = x(goodX);
    y = y(goodX);

    if numel(x) < 2
        error('plot_line1d:TooFewPoints', 'At least 2 valid x points are required.');
    end
end

function y = applyNaNMode(y, x, maskNaNRaw, nanMode)
    mode = lower(strtrim(nanMode));
    switch mode
        case 'break'
            y(maskNaNRaw) = NaN;
        case 'omit'
            y(maskNaNRaw) = NaN;
        case 'interp'
            m = ~maskNaNRaw & isfinite(y);
            if sum(m) >= 2
                y(maskNaNRaw) = interp1(x(m), y(m), x(maskNaNRaw), 'linear', 'extrap');
            else
                y(maskNaNRaw) = NaN;
            end
        otherwise
            error('plot_line1d:BadNaNMode', 'NaNMode must be break/omit/interp.');
    end
end

function ys = gaussianSmoothNaN1d(y, sigma)
    if sigma <= 0
        ys = y;
        return;
    end

    r = max(1, ceil(3*sigma));
    t = -r:r;
    g = exp(-(t.^2)/(2*sigma^2));
    g = g / sum(g);

    w = double(isfinite(y));
    y0 = y;
    y0(~isfinite(y0)) = 0;

    num = conv(y0, g, 'same');
    den = conv(w, g, 'same');

    ys = num ./ den;
    ys(den <= eps) = NaN;
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
