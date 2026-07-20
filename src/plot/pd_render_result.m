function handles = pd_render_result(ax, result, options, viewName)
%RENDERRESULT Render a POST_DATA numeric result into caller-owned axes.
% The function never creates a figure and is MATLAB R2016b compatible.

    if nargin < 2 || ~ishghandle(ax)
        error('postdata:renderResult:BadAxes', 'A valid target axes is required.');
    end
    if ~isstruct(result) || ~isfield(result, 'analysisType')
        error('postdata:renderResult:BadResult', ...
            'Result must contain analysisType metadata from postdata_run.');
    end

    if nargin < 3 || isempty(options)
        catalog = pd_plot_option_catalog();
        options = pd_plot_options_from_table(catalog, pd_catalog_table_data(catalog));
    end
    views = pd_result_plot_views(result);
    if nargin < 4 || isempty(viewName)
        viewName = views(1).Id;
    end
    viewName = lower(strtrim(pd_to_char(viewName)));
    if ~any(strcmp(viewName, {views.Id}))
        error('postdata:renderResult:BadView', ...
            'View "%s" is unavailable for this result.', viewName);
    end

    if ~isfield(options, 'UpdateMode')
        options.UpdateMode = 'replace';
    end
    updateMode = strtrim(pd_to_char(options.UpdateMode));
    overlay = strcmpi(updateMode, 'overlay');
    if ~overlay
        cla(ax, 'reset');
    end
    wasHeld = ishold(ax);
    if overlay
        hold(ax, 'on');
    end
    applyPrePlotStyle(ax, options);
    holdCleanup = onCleanup(@() restoreHoldState(ax, wasHeld, overlay));
    analysisType = lower(strtrim(pd_to_char(result.analysisType)));
    plotData = pd_result_plot_data(result, viewName);

    switch analysisType
        case 'chunk'
            handles = renderChunk(ax, result, options, viewName, plotData);
        case 'cluster'
            handles = renderCluster(ax, result, options, viewName, plotData);
        case 'vx'
            handles = renderVelocity(ax, result, options, viewName, plotData);
        case 'massx'
            handles = renderMassX(ax, result, options, viewName, plotData);
        case 'network2d'
            handles = renderNetwork(ax, result, options, viewName, plotData);
        otherwise
            error('postdata:renderResult:UnknownType', ...
                'Unsupported analysis type: %s', analysisType);
    end
    pd_apply_publication_style(ax, options);
end

function handles = renderMassX(ax, result, options, viewName, plotData)
    if isempty(result.cumulativeDensity)
        text(0.5, 0.5, 'No density variables', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        handles = [];
        return;
    end
    if strcmp(viewName, 'particle-count')
        titleText = 'SPH particle count';
        yLabel = 'Particle count';
    elseif strcmp(viewName, 'differential')
        titleText = 'SPH mass-x local';
        yLabel = 'Areal density (mg/cm^2)';
    else
        titleText = 'SPH mass-x cumulative';
        yLabel = 'Cumulative areal density (mg/cm^2)';
    end
    handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2:end), ...
        'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
        'Marker', options.MarkerSymbol, 'MarkerSize', options.MarkerSize);
    baseLabel = resultSeriesLabel(result);
    for i = 1:numel(handles)
        if i <= numel(result.densityVars)
            if strcmpi(options.UpdateMode, 'overlay')
                labelText = sprintf('%s: %s', baseLabel, result.densityVars{i});
            else
                labelText = result.densityVars{i};
            end
            set(handles(i), 'DisplayName', labelText);
        end
    end
    label = result.coordinateLabel;
    if isfield(result, 'coordinateUnit') && ~isempty(result.coordinateUnit)
        label = sprintf('%s (%s)', label, result.coordinateUnit);
    end
    xlabel(ax, label, 'Interpreter', 'none');
    ylabel(ax, yLabel);
    title(ax, sprintf('%s @ %g', titleText, result.timestep));
end

function handles = renderChunk(ax, result, options, viewName, plotData)
    if strcmp(viewName, 'histogram')
        handles = bar(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            1.0, 'FaceColor', firstPaletteColor(options));
        xlabel(ax, result.variableUsed, 'Interpreter', 'none');
        ylabel(ax, 'Count');
        title(ax, sprintf('%s histogram @ timestep %g', result.variableUsed, result.timestep), ...
            'Interpreter', 'none');
        return;
    elseif any(strcmp(viewName, {'profile-x','profile-y'}))
        if strcmp(viewName, 'profile-x'), label = 'x'; else, label = 'y'; end
        handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            'LineWidth', options.LineWidth, ...
            'LineStyle', options.LineStyle, 'Marker', options.MarkerSymbol, ...
            'MarkerSize', options.MarkerSize, ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, label);
        ylabel(ax, ['Mean ', result.variableUsed], 'Interpreter', 'none');
        title(ax, sprintf('%s mean profile %s @ timestep %g', ...
            result.variableUsed, upper(label), result.timestep), 'Interpreter', 'none');
        return;
    elseif isfield(result, 'dimension') && strcmpi(result.dimension, '3d')
        handles = scatter3(ax, plotData.Values(:, 1), ...
            plotData.Values(:, 2), plotData.Values(:, 3), ...
            options.MarkerSize ^ 2, plotData.Values(:, 4), ...
            scatterMarker(options.MarkerSymbol), 'filled', ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        zlabel(ax, 'z');
        axis(ax, 'tight');
        grid(ax, 'on');
        setColorbar(ax, options.ShowColorbar, result.variableUsed);
    elseif isempty(result.y)
        handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            'LineWidth', options.LineWidth, ...
            'LineStyle', options.LineStyle, 'Marker', options.MarkerSymbol, ...
            'MarkerSize', options.MarkerSize, ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, 'x');
        ylabel(ax, result.variableUsed, 'Interpreter', 'none');
    else
        handles = renderField2D(ax, result, options, plotData);
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        axis(ax, 'tight');
        setColorbar(ax, options.ShowColorbar, result.variableUsed);
    end
    title(ax, sprintf('%s @ timestep %g', result.variableUsed, result.timestep), ...
        'Interpreter', 'none');
end

function handles = renderField2D(ax, result, options, plotData)
    isGrid = plotData.IsRectangularGrid;
    xValues = plotData.GridX;
    yValues = plotData.GridY;
    valueGrid = plotData.GridZ;
    mode = lower(strtrim(pd_to_char(options.FieldRenderMode)));
    if strcmp(mode, 'auto')
        if isGrid, mode = 'image'; else, mode = 'scatter'; end
    end
    if any(strcmp(mode, {'image','contour'})) && ~isGrid
        if strcmp(mode, 'image')
            warning('postdata:IrregularFieldFallback', ...
                'FieldRenderMode=image requires a rectangular grid; using scatter.');
            mode = 'scatter';
        else
            error('postdata:ContourRequiresGrid', ...
                'FieldRenderMode=contour requires a complete rectangular grid.');
        end
    end
    switch mode
        case 'image'
            handles = imagesc(xValues, yValues, valueGrid, 'Parent', ax);
            set(ax, 'YDir', 'normal');
            if strcmpi(options.UpdateMode, 'overlay')
                set(handles, 'AlphaData', 0.48 .* double(isfinite(valueGrid)));
            end
            addLegendProxy(ax, resultSeriesLabel(result));
        case 'contour'
            [~, handles] = contourf(ax, xValues, yValues, valueGrid, ...
                round(options.ContourLevels), 'LineStyle', 'none');
            if strcmpi(options.UpdateMode, 'overlay')
                try
                    set(handles, 'FaceAlpha', 0.48);
                catch
                    % FaceAlpha support differs between graphics releases.
                end
            end
            addLegendProxy(ax, resultSeriesLabel(result));
        otherwise
            handles = scatter(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
                options.MarkerSize ^ 2, plotData.Values(:, 3), ...
                scatterMarker(options.MarkerSymbol), 'filled', ...
                'DisplayName', resultSeriesLabel(result));
    end
end

function handle = addLegendProxy(ax, label)
    wasHeld = ishold(ax);
    hold(ax, 'on');
    handle = plot(ax, NaN, NaN, 'DisplayName', label);
    if ~wasHeld, hold(ax, 'off'); end
end

function handles = renderCluster(ax, result, options, viewName, plotData)
    if isempty(result.diameter)
        text(0.5, 0.5, 'No valid clusters', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        handles = [];
        return;
    end
    if strcmp(viewName, 'cdf')
        handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
            'Marker', options.MarkerSymbol, 'MarkerSize', options.MarkerSize);
        xlabel(ax, 'Equivalent diameter'); ylabel(ax, 'Cumulative probability');
        title(ax, sprintf('Cluster diameter CDF @ timestep %g', result.timestep));
        return;
    elseif strcmp(viewName, 'mean')
        handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
            'Marker', effectiveMarker(options.MarkerSymbol, 'o'), ...
            'MarkerSize', options.MarkerSize);
        xlabel(ax, 'Position bin'); ylabel(ax, 'Mean equivalent diameter');
        title(ax, sprintf('Mean cluster diameter @ timestep %g', result.timestep));
        return;
    elseif strcmp(viewName, 'probability')
        yLabel = 'Probability';
    else
        yLabel = 'Count';
    end
    handles = bar(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
        1.0, 'FaceColor', firstPaletteColor(options));
    set(handles, 'LineWidth', options.LineWidth);
    set(handles, 'DisplayName', resultSeriesLabel(result));
    if strcmpi(options.UpdateMode, 'overlay')
        try
            set(handles, 'FaceAlpha', 0.35);
        catch
            % FaceAlpha is renderer-dependent on older MATLAB releases.
        end
    end
    xlabel(ax, 'Equivalent diameter');
    ylabel(ax, yLabel);
    title(ax, sprintf('Cluster distribution @ timestep %g', result.timestep));
end

function handles = renderVelocity(ax, result, options, viewName, plotData)
    if isempty(result.cumulativeDensity)
        text(0.5, 0.5, 'No density variables', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        handles = [];
        return;
    end
    if strcmp(viewName, 'differential')
        yLabel = 'Areal density (mg/cm^2)';
        titleText = 'mass-v differential distribution';
    else
        yLabel = 'Cumulative areal density (mg/cm^2)';
        titleText = 'mass-v cumulative distribution';
    end
    handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2:end), ...
        'LineWidth', options.LineWidth, ...
        'LineStyle', options.LineStyle, 'Marker', options.MarkerSymbol, ...
        'MarkerSize', options.MarkerSize);
    baseLabel = resultSeriesLabel(result);
    for i = 1:numel(handles)
        if i <= numel(result.densityVars)
            if strcmpi(options.UpdateMode, 'overlay')
                labelText = sprintf('%s: %s', baseLabel, result.densityVars{i});
            else
                labelText = result.densityVars{i};
            end
            set(handles(i), 'DisplayName', labelText);
        else
            set(handles(i), 'DisplayName', baseLabel);
        end
    end
    label = 'Velocity';
    if isfield(result, 'velocityLabel'), label = result.velocityLabel; end
    if isfield(result, 'velocityUnit') && ~isempty(result.velocityUnit)
        label = sprintf('%s (%s)', label, result.velocityUnit);
    end
    xlabel(ax, label, 'Interpreter', 'none');
    ylabel(ax, yLabel);
    title(ax, sprintf('%s @ timestep %g', titleText, result.timestep));
end

function handles = renderNetwork(ax, result, options, viewName, plotData)
    if any(strcmp(viewName, {'pore-label','matrix-label'}))
        if strcmp(viewName, 'pore-label'), phaseName = 'pore'; else, phaseName = 'matrix'; end
        handles = imagesc(plotData.GridX, plotData.GridY, ...
            plotData.GridZ, 'Parent', ax);
        set(ax, 'YDir', 'normal'); axis(ax, 'tight');
        setColorbar(ax, options.ShowColorbar, [phaseName, ' component label']);
        xlabel(ax, 'x'); ylabel(ax, 'y');
        title(ax, sprintf('%s components @ timestep %g', phaseName, result.timestep));
        return;
    elseif any(strcmp(viewName, {'pore-diameter','matrix-diameter'}))
        if strcmp(viewName, 'pore-diameter'), phaseName = 'pore'; else, phaseName = 'matrix'; end
        handles = bar(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            1.0, 'FaceColor', firstPaletteColor(options));
        xlabel(ax, 'Equivalent diameter'); ylabel(ax, 'Count');
        title(ax, sprintf('%s diameter @ timestep %g', phaseName, result.timestep));
        return;
    elseif any(strcmp(viewName, {'profile-x','profile-y'}))
        axisName = viewName(end);
        handles = plot(ax, plotData.Values(:, 1), plotData.Values(:, 2), ...
            'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
            'Marker', options.MarkerSymbol, 'MarkerSize', options.MarkerSize, ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, axisName); ylabel(ax, 'Porosity');
        title(ax, sprintf('Porosity profile %s @ timestep %g', ...
            upper(axisName), result.timestep));
        return;
    end
    phase = plotData.GridZ;
    handles = imagesc(plotData.GridX, plotData.GridY, phase, 'Parent', ax);
    if strcmpi(options.UpdateMode, 'overlay')
        set(handles, 'AlphaData', 0.45 .* double(isfinite(phase)));
        plot(ax, NaN, NaN, 'LineWidth', options.LineWidth, ...
            'DisplayName', resultSeriesLabel(result));
    end
    set(ax, 'YDir', 'normal');
    axis(ax, 'tight');
    set(ax, 'CLim', [0, 1]);
    setColorbar(ax, options.ShowColorbar, 'Pore fraction');
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    title(ax, sprintf('Pore fraction @ timestep %g', result.timestep));
end

function applyPrePlotStyle(ax, options)
    set(ax, 'ColorOrder', pd_resolve_plot_palette(options));
end

function color = firstPaletteColor(options)
    colors = pd_resolve_plot_palette(options);
    color = colors(1, :);
end

function marker = scatterMarker(value)
    marker = pd_to_char(value);
    if strcmpi(marker, 'none'), marker = 'o'; end
end

function marker = effectiveMarker(value, fallback)
    marker = pd_to_char(value);
    if strcmpi(marker, 'none'), marker = fallback; end
end

function setColorbar(ax, enabled, label)
    if nargin < 3, label = ''; end
    if enabled
        cb = colorbar('peer', ax);
        if ~isempty(label), ylabel(cb, label, 'Interpreter', 'none'); end
    else
        colorbar('peer', ax, 'off');
    end
end

function label = resultSeriesLabel(result)
    base = lower(strtrim(pd_to_char(result.analysisType)));
    if isfield(result, 'filePath') && ~isempty(result.filePath)
        [~, name, extension] = fileparts(pd_to_char(result.filePath));
        base = [name, extension];
    end
    if isfield(result, 'timestep') && isscalar(result.timestep) && isfinite(result.timestep)
        label = sprintf('%s @ %g', base, result.timestep);
    else
        label = base;
    end
end

function restoreHoldState(ax, wasHeld, overlay)
    if ~overlay || ~ishghandle(ax)
        return;
    end
    if wasHeld
        hold(ax, 'on');
    else
        hold(ax, 'off');
    end
end
