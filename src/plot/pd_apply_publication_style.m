function options = pd_apply_publication_style(ax, options)
%PD_APPLY_PUBLICATION_STYLE Apply the shared editable paper-figure style.
%   OPTIONS may be the complete plot-options structure, a partial structure,
%   or a name/value cell array. Existing titles and labels are preserved when
%   the corresponding override option is empty.

    if nargin < 1 || isempty(ax) || ~all(ishghandle(ax))
        error('postdata:BadStyleAxes', 'Valid axes are required.');
    end
    if nargin < 2, options = []; end
    options = pd_normalize_plot_options(options);
    for k = 1:numel(ax)
        styleOneAxes(ax(k), options);
    end
end

function styleOneAxes(ax, options)
    palette = pd_resolve_plot_palette(options);
    set(ax, 'ColorOrder', palette, 'FontName', options.FontName, ...
        'FontSize', options.FontSize, 'LineWidth', options.AxisLineWidth, ...
        'TickDir', options.TickDirection, 'GridLineStyle', ...
        options.GridLineStyle, 'Layer', 'top');
    try, set(ax, 'GridAlpha', options.GridAlpha); catch, end
    if options.ShowGrid, grid(ax, 'on'); else, grid(ax, 'off'); end
    if options.ShowMinorGrid
        set(ax, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    else
        set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    end
    if options.BoxOn, box(ax, 'on'); else, box(ax, 'off'); end
    if options.AxisEqual, axis(ax, 'equal'); end
    if options.LogX, set(ax, 'XScale', 'log'); else, set(ax, 'XScale', 'linear'); end
    if options.LogY, set(ax, 'YScale', 'log'); else, set(ax, 'YScale', 'linear'); end
    if options.ReverseX, set(ax, 'XDir', 'reverse'); else, set(ax, 'XDir', 'normal'); end
    if options.ReverseY, set(ax, 'YDir', 'reverse'); else, set(ax, 'YDir', 'normal'); end

    applyRange(ax, 'XLim', options.XLim);
    applyRange(ax, 'YLim', options.YLim);
    applyRange(ax, 'CLim', options.ColorLimits);
    applyColormap(ax, options.Colormap);
    styleSeries(ax, options, palette);
    updateColorbar(ax, options);

    if ~isempty(options.Title), title(ax, options.Title); end
    if ~isempty(options.XLabel), xlabel(ax, options.XLabel); end
    if ~isempty(options.YLabel), ylabel(ax, options.YLabel); end
    set(get(ax, 'Title'), 'FontName', options.FontName, ...
        'FontSize', options.TitleFontSize, 'FontWeight', options.TitleWeight, ...
        'Interpreter', options.TextInterpreter);
    set(get(ax, 'XLabel'), 'FontName', options.FontName, ...
        'FontSize', options.FontSize, 'Interpreter', options.TextInterpreter);
    set(get(ax, 'YLabel'), 'FontName', options.FontName, ...
        'FontSize', options.FontSize, 'Interpreter', options.TextInterpreter);
    applyLegendStyle(ax, options);
end

function updateColorbar(ax, options)
    colored = ~isempty(findobj(ax, 'Type', 'image')) || ...
        ~isempty(findobj(ax, 'Type', 'surface')) || ...
        ~isempty(findobj(ax, 'Type', 'scatter')) || ...
        ~isempty(findobj(ax, 'Type', 'patch'));
    if options.ShowColorbar && colored
        cb = colorbar('peer', ax);
        try
            set(cb, 'FontName', options.FontName, 'FontSize', options.FontSize);
        catch
            % Colorbar font properties differ between graphics releases.
        end
        if ~isempty(options.ColorbarLabel)
            ylabel(cb, options.ColorbarLabel, 'Interpreter', options.TextInterpreter);
        end
    elseif ~options.ShowColorbar
        colorbar('peer', ax, 'off');
    end
end

function styleSeries(ax, options, palette)
    linesFound = findobj(ax, 'Type', 'line');
    linesFound = linesFound(end:-1:1);
    for i = 1:numel(linesFound)
        [color, lineStyle] = seriesAppearance(i, options, palette);
        set(linesFound(i), 'LineWidth', options.LineWidth, ...
            'LineStyle', lineStyle, 'Marker', options.MarkerSymbol, ...
            'MarkerSize', options.MarkerSize, 'Color', color);
    end
    bars = findobj(ax, 'Type', 'bar');
    for i = 1:numel(bars)
        try
            set(bars(i), 'FaceColor', palette(mod(i - 1, size(palette, 1)) + 1, :), ...
                'LineWidth', options.AxisLineWidth);
        catch
            % Bar property support differs between HG releases.
        end
    end
end

function [color, lineStyle] = seriesAppearance(index, options, palette)
    styles = lineStyleCycle(options.LineStyle);
    lineStyle = options.LineStyle;
    switch lower(pd_to_char(options.SeriesStyleMode))
        case 'monochrome'
            color = [0.08 0.08 0.08];
            lineStyle = styles{mod(index - 1, numel(styles)) + 1};
        case 'color-and-style'
            color = palette(mod(index - 1, size(palette, 1)) + 1, :);
            lineStyle = styles{mod(index - 1, numel(styles)) + 1};
        otherwise
            color = palette(mod(index - 1, size(palette, 1)) + 1, :);
    end
end

function styles = lineStyleCycle(firstStyle)
    candidates = {firstStyle, '-', '--', ':', '-.'};
    styles = {};
    for i = 1:numel(candidates)
        if ~any(strcmp(candidates{i}, styles))
            styles{end + 1} = candidates{i}; %#ok<AGROW>
        end
    end
end

function applyLegendStyle(ax, options)
    if ~options.ShowLegend
        legend(ax, 'off');
        return;
    end
    objects = findobj(ax, '-property', 'DisplayName');
    namedCount = 0;
    for i = 1:numel(objects)
        name = get(objects(i), 'DisplayName');
        if ischar(name) && ~isempty(name)
            namedCount = namedCount + 1;
        end
    end
    if namedCount < 2
        legend(ax, 'off');
        return;
    end
    lgd = legend(ax, 'show', 'Interpreter', options.TextInterpreter, ...
        'Location', options.LegendLocation);
    if options.LegendBox, set(lgd, 'Box', 'on'); else, set(lgd, 'Box', 'off'); end
    set(lgd, 'FontName', options.FontName, ...
        'FontSize', max(8, options.FontSize - 1));
end

function applyRange(ax, propertyName, value)
    if isempty(value), return; end
    if ~(isnumeric(value) && numel(value) == 2 && all(isfinite(value)) && ...
            value(2) > value(1))
        error('postdata:BadPlotRange', ...
            '%s must be empty or [min max].', propertyName);
    end
    set(ax, propertyName, value(:).');
end

function applyColormap(ax, name)
    valid = {'parula','jet','hsv','hot','cool','spring','summer','autumn', ...
        'winter','gray','bone','copper','pink','lines','colorcube','prism', ...
        'flag','white'};
    name = lower(strtrim(pd_to_char(name)));
    if ~any(strcmp(name, valid))
        error('postdata:BadColormap', 'Unsupported colormap: %s.', name);
    end
    mapFunction = str2func(name);
    colormap(ax, mapFunction(256));
end
