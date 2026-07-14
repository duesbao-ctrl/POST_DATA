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

    switch analysisType
        case 'chunk'
            handles = renderChunk(ax, result, options, viewName);
        case 'cluster'
            handles = renderCluster(ax, result, options, viewName);
        case 'vx'
            handles = renderVelocity(ax, result, options, viewName);
        case 'massx'
            handles = renderMassX(ax, result, options, viewName);
        case 'network2d'
            handles = renderNetwork(ax, result, options, viewName);
        otherwise
            error('postdata:renderResult:UnknownType', ...
                'Unsupported analysis type: %s', analysisType);
    end
    pd_apply_publication_style(ax, options);
end

function handles = renderMassX(ax, result, options, viewName)
    if isempty(result.cumulativeDensity)
        text(0.5, 0.5, 'No density variables', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        handles = [];
        return;
    end
    if strcmp(viewName, 'differential')
        values = result.density;
        titleText = 'mass-x differential distribution';
        yLabel = 'Areal density (mg/cm^2)';
    else
        values = result.cumulativeDensity;
        titleText = 'mass-x cumulative distribution';
        yLabel = 'Cumulative areal density (mg/cm^2)';
    end
    handles = plot(ax, result.coordinate, values, ...
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
    title(ax, sprintf('%s @ timestep %g', titleText, result.timestep));
end

function handles = renderChunk(ax, result, options, viewName)
    if strcmp(viewName, 'histogram')
        values = result.value(isfinite(result.value));
        [counts, centers] = hist(values, max(5, min(50, round(sqrt(numel(values)))))); %#ok<HIST>
        handles = bar(ax, centers, counts, 1.0, 'FaceColor', firstPaletteColor(options));
        xlabel(ax, result.variableUsed, 'Interpreter', 'none');
        ylabel(ax, 'Count');
        title(ax, sprintf('%s histogram @ timestep %g', result.variableUsed, result.timestep), ...
            'Interpreter', 'none');
        return;
    elseif any(strcmp(viewName, {'profile-x','profile-y'}))
        if strcmp(viewName, 'profile-x'), coordinate = result.x; label = 'x'; ...
        else, coordinate = result.y; label = 'y'; end
        [centers, means] = coordinateMean(coordinate, result.value);
        handles = plot(ax, centers, means, 'LineWidth', options.LineWidth, ...
            'LineStyle', options.LineStyle, 'Marker', options.MarkerSymbol, ...
            'MarkerSize', options.MarkerSize, ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, label);
        ylabel(ax, ['Mean ', result.variableUsed], 'Interpreter', 'none');
        title(ax, sprintf('%s mean profile %s @ timestep %g', ...
            result.variableUsed, upper(label), result.timestep), 'Interpreter', 'none');
        return;
    elseif isempty(result.y)
        handles = plot(ax, result.x, result.value, 'LineWidth', options.LineWidth, ...
            'LineStyle', options.LineStyle, 'Marker', options.MarkerSymbol, ...
            'MarkerSize', options.MarkerSize, ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, 'x');
        ylabel(ax, result.variableUsed, 'Interpreter', 'none');
    else
        handles = renderField2D(ax, result, options);
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        axis(ax, 'tight');
        setColorbar(ax, options.ShowColorbar, result.variableUsed);
    end
    title(ax, sprintf('%s @ timestep %g', result.variableUsed, result.timestep), ...
        'Interpreter', 'none');
end

function handles = renderField2D(ax, result, options)
    [isGrid, xValues, yValues, valueGrid] = rectangularGrid( ...
        result.x, result.y, result.value);
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
                try, set(handles, 'FaceAlpha', 0.48); catch, end
            end
            addLegendProxy(ax, resultSeriesLabel(result));
        otherwise
            handles = scatter(ax, result.x, result.y, options.MarkerSize ^ 2, ...
                result.value, scatterMarker(options.MarkerSymbol), 'filled', ...
                'DisplayName', resultSeriesLabel(result));
    end
end

function handle = addLegendProxy(ax, label)
    wasHeld = ishold(ax);
    hold(ax, 'on');
    handle = plot(ax, NaN, NaN, 'DisplayName', label);
    if ~wasHeld, hold(ax, 'off'); end
end

function [isGrid, xValues, yValues, valueGrid] = rectangularGrid(x, y, value)
    x = x(:); y = y(:); value = value(:);
    [xValues, ~, xIndex] = unique(x);
    [yValues, ~, yIndex] = unique(y);
    isGrid = numel(xValues) * numel(yValues) == numel(value);
    valueGrid = [];
    if ~isGrid, return; end
    linearIndex = sub2ind([numel(yValues), numel(xValues)], yIndex, xIndex);
    if numel(unique(linearIndex)) ~= numel(linearIndex)
        isGrid = false;
        return;
    end
    valueGrid = nan(numel(yValues), numel(xValues));
    valueGrid(linearIndex) = value;
end

function [centers, means] = coordinateMean(coordinate, values)
    valid = isfinite(coordinate) & isfinite(values);
    [centers, ~, groups] = unique(coordinate(valid));
    means = accumarray(groups, values(valid), [], @mean);
end

function handles = renderCluster(ax, result, options, viewName)
    if isempty(result.diameter)
        text(0.5, 0.5, 'No valid clusters', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        handles = [];
        return;
    end
    if strcmp(viewName, 'cdf')
        values = sort(result.diameter(:));
        handles = plot(ax, values, (1:numel(values)).' ./ numel(values), ...
            'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
            'Marker', options.MarkerSymbol, 'MarkerSize', options.MarkerSize);
        xlabel(ax, 'Equivalent diameter'); ylabel(ax, 'Cumulative probability');
        title(ax, sprintf('Cluster diameter CDF @ timestep %g', result.timestep));
        return;
    elseif strcmp(viewName, 'mean')
        handles = plot(ax, result.meanByBin.centers, result.meanByBin.meanDiameter, ...
            'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
            'Marker', effectiveMarker(options.MarkerSymbol, 'o'), ...
            'MarkerSize', options.MarkerSize);
        xlabel(ax, 'Position bin'); ylabel(ax, 'Mean equivalent diameter');
        title(ax, sprintf('Mean cluster diameter @ timestep %g', result.timestep));
        return;
    elseif strcmp(viewName, 'probability') && isfield(result, 'hist')
        centers = result.hist.centers;
        counts = result.hist.prob;
        yLabel = 'Probability';
    else
        nBins = max(5, min(50, round(sqrt(numel(result.diameter)))));
        [counts, centers] = hist(result.diameter, nBins); %#ok<HIST>
        yLabel = 'Count';
    end
    handles = bar(ax, centers, counts, 1.0, 'FaceColor', firstPaletteColor(options));
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

function handles = renderVelocity(ax, result, options, viewName)
    if isempty(result.cumulativeDensity)
        text(0.5, 0.5, 'No density variables', 'Parent', ax, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        handles = [];
        return;
    end
    if strcmp(viewName, 'differential')
        values = result.density;
        yLabel = 'Areal density (mg/cm^2)';
        titleText = 'mass-v differential distribution';
    else
        values = result.cumulativeDensity;
        yLabel = 'Cumulative areal density (mg/cm^2)';
        titleText = 'mass-v cumulative distribution';
    end
    handles = plot(ax, result.velocity, values, 'LineWidth', options.LineWidth, ...
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

function handles = renderNetwork(ax, result, options, viewName)
    if any(strcmp(viewName, {'pore-label','matrix-label'}))
        if strcmp(viewName, 'pore-label'), phaseName = 'pore'; else, phaseName = 'matrix'; end
        handles = imagesc(result.xCenters, result.yCenters, ...
            result.(phaseName).labelGrid, 'Parent', ax);
        set(ax, 'YDir', 'normal'); axis(ax, 'tight');
        setColorbar(ax, options.ShowColorbar, [phaseName, ' component label']);
        xlabel(ax, 'x'); ylabel(ax, 'y');
        title(ax, sprintf('%s components @ timestep %g', phaseName, result.timestep));
        return;
    elseif any(strcmp(viewName, {'pore-diameter','matrix-diameter'}))
        if strcmp(viewName, 'pore-diameter'), phaseName = 'pore'; else, phaseName = 'matrix'; end
        values = result.(phaseName).components.equivDiameter;
        [counts, centers] = hist(values, max(3, min(40, round(sqrt(numel(values)))))); %#ok<HIST>
        handles = bar(ax, centers, counts, 1.0, 'FaceColor', firstPaletteColor(options));
        xlabel(ax, 'Equivalent diameter'); ylabel(ax, 'Count');
        title(ax, sprintf('%s diameter @ timestep %g', phaseName, result.timestep));
        return;
    elseif any(strcmp(viewName, {'profile-x','profile-y'}))
        axisName = viewName(end);
        profile = result.profile.(axisName);
        handles = plot(ax, profile.centers, profile.porosity, ...
            'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
            'Marker', options.MarkerSymbol, 'MarkerSize', options.MarkerSize, ...
            'DisplayName', resultSeriesLabel(result));
        xlabel(ax, axisName); ylabel(ax, 'Porosity');
        title(ax, sprintf('Porosity profile %s @ timestep %g', ...
            upper(axisName), result.timestep));
        return;
    end
    phase = double(result.poreMask);
    phase(~result.validMask) = NaN;
    if isfield(result, 'cutCell') && isfield(result.cutCell, 'pore') && ...
            isfield(result.cutCell.pore, 'fraction')
        phase = result.cutCell.pore.fraction;
    end
    handles = imagesc(result.xCenters, result.yCenters, phase, 'Parent', ax);
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
