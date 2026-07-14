function plots = pd_render_network2d_snapshot(xCenters, yCenters, phaseGrid, ...
        validMask, pore, matrix, profile, cutCell, timestep, opt)
%PD_RENDER_NETWORK2D_SNAPSHOT Render all non-evolution network2d figures.

    plots = emptyPlotsStruct();
    if ~opt.MakePlots
        return;
    end

    plots.phaseFig = plotPhaseGrid(xCenters, yCenters, phaseGrid, timestep, ...
        opt.ThresholdN, opt.PlotRangeX, opt.PlotRangeY, opt.CutCellPlotRefinement);
    plots.poreLabelFig = plotLabelGrid(xCenters, yCenters, ...
        pore.components.ownerGrid, validMask, pore.sizeRank, ...
        sprintf('Pore Labels @ timestep %g', timestep), ...
        opt.PlotRangeX, opt.PlotRangeY, cutCell.pore.fraction, ...
        opt.CutCellPlotRefinement);
    plots.matrixLabelFig = plotLabelGrid(xCenters, yCenters, ...
        matrix.components.ownerGrid, validMask, matrix.sizeRank, ...
        sprintf('Matrix Labels @ timestep %g', timestep), ...
        opt.PlotRangeX, opt.PlotRangeY, cutCell.matrix.fraction, ...
        opt.CutCellPlotRefinement);
    plots.connectivityFig = plotConnectivityHighlights( ...
        xCenters, yCenters, validMask, pore, matrix, timestep, ...
        opt.PlotRangeX, opt.PlotRangeY);
    plots.poreCountFig = plotCountDistribution( ...
        pore.equivDiameterDistribution, 'Pore', opt.DiameterRange, ...
        opt.HistScale, opt.DiameterPlotStyle);
    plots.matrixCountFig = plotCountDistribution( ...
        matrix.equivDiameterDistribution, 'Matrix', opt.DiameterRange, ...
        opt.HistScale, opt.DiameterPlotStyle);

    if strcmp(opt.PositionAxis, 'x') || strcmp(opt.PositionAxis, 'both')
        plots.porePositionFigX = plotPositionDistribution( ...
            pore.positionDistribution.x, 'Pore', 'x');
        plots.matrixPositionFigX = plotPositionDistribution( ...
            matrix.positionDistribution.x, 'Matrix', 'x');
    end
    if strcmp(opt.PositionAxis, 'y') || strcmp(opt.PositionAxis, 'both')
        plots.porePositionFigY = plotPositionDistribution( ...
            pore.positionDistribution.y, 'Pore', 'y');
        plots.matrixPositionFigY = plotPositionDistribution( ...
            matrix.positionDistribution.y, 'Matrix', 'y');
    end
    if isfield(profile, 'x')
        plots.profileFigX = plotDirectionalProfile(profile.x, timestep);
    end
    if isfield(profile, 'y')
        plots.profileFigY = plotDirectionalProfile(profile.y, timestep);
    end
end
function fig = plotPhaseGrid(xCenters, yCenters, phaseGrid, timestep, thresholdN, xRange, yRange, refinement)
    if nargin < 8 || isempty(refinement)
        refinement = 1;
    end
    fig = figure('Color', 'w', 'Name', '2D Network Phase Map');
    ax = axes('Parent', fig);
    [plotX, plotY, plotZ, plotAlpha] = refinePhaseGridForPlot(xCenters, yCenters, phaseGrid, refinement);
    hImg = imagesc(ax, plotX, plotY, plotZ);
    set(hImg, 'AlphaData', plotAlpha);
    set(ax, 'Color', [0.72 0.72 0.72]);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    applyPlotRanges(ax, xCenters, yCenters, xRange, yRange);
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
    colormap(ax, buildPhaseFractionColormap());
    setAxesCLim(ax, [0 1]);
    cb = colorbar(ax, 'Ticks', [0, 0.5, 1], ...
        'TickLabels', {'Matrix', 'Mixed', 'Pore'});
    ylabel(cb, 'Pore area fraction');
    xlabel(ax, 'x', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'y', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, sprintf('Phase map (Ncount < %g is pore) @ timestep %g', thresholdN, timestep), ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = plotLabelGrid(xCenters, yCenters, labelGrid, validMask, sizeRank, ttl, xRange, yRange, phaseFraction, refinement)
    if nargin < 9
        phaseFraction = [];
    end
    if nargin < 10 || isempty(refinement)
        refinement = 1;
    end
    fig = figure('Color', 'w', 'Name', ttl);
    ax = axes('Parent', fig);
    z = buildSizeRankLabelGrid(labelGrid, validMask, sizeRank);
    [plotX, plotY, plotZ, plotAlpha] = refineLabelGridForPlot( ...
        xCenters, yCenters, z, validMask, phaseFraction, refinement);
    hImg = imagesc(ax, plotX, plotY, plotZ);
    set(hImg, 'AlphaData', plotAlpha);
    set(ax, 'Color', [0.93 0.93 0.91]);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    applyPlotRanges(ax, xCenters, yCenters, xRange, yRange);
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
    cmap = buildSizeRankColormap(numel(sizeRank.rank));
    colormap(ax, cmap);
    setAxesCLim(ax, [0 max(numel(sizeRank.rank), 1)]);
    cb = colorbar(ax);
    ylabel(cb, 'Size Rank (largest first)');
    xlabel(ax, 'x', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'y', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, ttl, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    set(fig, 'UserData', sizeRank);
end

function [plotX, plotY, plotZ, plotAlpha] = refineLabelGridForPlot( ...
        xCenters, yCenters, rankGrid, validMask, phaseFraction, refinement)
    refinement = max(1, round(refinement));
    if isempty(phaseFraction)
        phaseFraction = double(rankGrid > 0);
        phaseFraction(~validMask) = NaN;
    end

    plotX = xCenters;
    plotY = yCenters;
    plotZ = rankGrid;
    plotZ(~isfinite(plotZ)) = 0;
    plotAlpha = phaseFraction;
    plotAlpha(~isfinite(plotAlpha)) = 0;
    plotAlpha = min(max(plotAlpha, 0), 1);
    plotAlpha = enhanceLabelAlphaForDisplay(plotAlpha);

    if refinement <= 1 || numel(xCenters) < 2 || numel(yCenters) < 2
        return;
    end

    nx = numel(xCenters);
    ny = numel(yCenters);
    plotX = linspace(xCenters(1), xCenters(end), (nx - 1) * refinement + 1);
    plotY = linspace(yCenters(1), yCenters(end), (ny - 1) * refinement + 1);
    [xGrid, yGrid] = meshgrid(xCenters, yCenters);
    [xQuery, yQuery] = meshgrid(plotX, plotY);

    rankForInterp = rankGrid;
    rankForInterp(~isfinite(rankForInterp)) = 0;
    plotZ = interp2(xGrid, yGrid, rankForInterp, xQuery, yQuery, 'nearest');
    plotZ(~isfinite(plotZ)) = 0;

    alphaForInterp = phaseFraction;
    alphaForInterp(~isfinite(alphaForInterp)) = 0;
    plotAlpha = interp2(xGrid, yGrid, alphaForInterp, xQuery, yQuery, 'linear');
    plotAlpha(~isfinite(plotAlpha)) = 0;
    plotAlpha = min(max(plotAlpha, 0), 1);
    plotAlpha = enhanceLabelAlphaForDisplay(plotAlpha);
end

function alpha = enhanceLabelAlphaForDisplay(alpha)
    visible = alpha > 0;
    alpha(visible) = 0.18 + 0.82 .* sqrt(alpha(visible));
    alpha(~visible) = 0;
end

function [plotX, plotY, plotZ, plotAlpha] = refinePhaseGridForPlot(xCenters, yCenters, phaseGrid, refinement)
    refinement = max(1, round(refinement));
    plotX = xCenters;
    plotY = yCenters;
    plotZ = phaseGrid;
    plotAlpha = double(isfinite(phaseGrid));
    plotZ(~isfinite(plotZ)) = 0;
    if refinement <= 1 || numel(xCenters) < 2 || numel(yCenters) < 2
        return;
    end

    nx = numel(xCenters);
    ny = numel(yCenters);
    plotX = linspace(xCenters(1), xCenters(end), (nx - 1) * refinement + 1);
    plotY = linspace(yCenters(1), yCenters(end), (ny - 1) * refinement + 1);
    [xGrid, yGrid] = meshgrid(xCenters, yCenters);
    [xQuery, yQuery] = meshgrid(plotX, plotY);
    valid = double(isfinite(phaseGrid));
    z = phaseGrid;
    z(~isfinite(z)) = 0;
    plotZ = interp2(xGrid, yGrid, z, xQuery, yQuery, 'linear');
    plotZ(~isfinite(plotZ)) = 0;
    plotZ = min(max(plotZ, 0), 1);
    plotAlpha = interp2(xGrid, yGrid, valid, xQuery, yQuery, 'linear');
    plotAlpha = double(plotAlpha > 0.999);
end

function cmap = buildPhaseFractionColormap()
    nPhase = 256;
    matrixColor = [0.18 0.38 0.62];
    mixedColor = [0.93 0.88 0.68];
    poreColor = [0.88 0.34 0.10];
    t = linspace(0, 1, nPhase).';
    colors = zeros(nPhase, 3);
    for i = 1:nPhase
        if t(i) <= 0.5
            a = t(i) / 0.5;
            colors(i, :) = (1 - a) * matrixColor + a * mixedColor;
        else
            a = (t(i) - 0.5) / 0.5;
            colors(i, :) = (1 - a) * mixedColor + a * poreColor;
        end
    end
    cmap = colors;
end

function z = buildSizeRankLabelGrid(labelGrid, validMask, sizeRank)
    z = zeros(size(labelGrid));
    z(~validMask) = NaN;
    if isempty(sizeRank.rankByLabel)
        return;
    end

    labels = labelGrid(validMask);
    mapped = zeros(size(labels));
    positive = labels > 0 & labels <= numel(sizeRank.rankByLabel);
    mapped(positive) = sizeRank.rankByLabel(labels(positive));
    z(validMask) = mapped;
end

function cmap = buildSizeRankColormap(numComponents)
    backgroundColor = [0.93 0.93 0.91];
    if numComponents <= 0
        cmap = [backgroundColor; backgroundColor];
        return;
    end

    base = [ ...
        0.88 0.32 0.10; ...
        0.19 0.58 0.33; ...
        0.78 0.53 0.12; ...
        0.70 0.21 0.58; ...
        0.12 0.61 0.66; ...
        0.56 0.40 0.72; ...
        0.72 0.25 0.26; ...
        0.37 0.60 0.18];

    if numComponents <= size(base, 1)
        colors = base(1:numComponents, :);
    else
        colors = [base; generateExtraComponentColors(numComponents - size(base, 1))];
    end
    cmap = [backgroundColor; colors];
end

function colors = generateExtraComponentColors(n)
    if n <= 0
        colors = zeros(0, 3);
        return;
    end

    k = (0:n-1).';
    hues = mod(0.04 + 0.61803398875 .* k, 1);
    blueBand = hues > 0.54 & hues < 0.70;
    hues(blueBand) = mod(hues(blueBand) + 0.20, 1);
    colors = hsv2rgb([hues, 0.68 * ones(n, 1), 0.86 * ones(n, 1)]);
end

function fig = plotConnectivityHighlights(xCenters, yCenters, validMask, pore, matrix, timestep, xRange, yRange)
    fig = figure('Color', 'w', 'Name', 'Directional Connectivity Highlights');

    subplot(2, 2, 1);
    drawConnectivityPanel(xCenters, yCenters, validMask, pore.mask, pore.labelGrid, ...
        pore.connectivity.x.componentIds, 'Pore connected in x', xRange, yRange);

    subplot(2, 2, 2);
    drawConnectivityPanel(xCenters, yCenters, validMask, pore.mask, pore.labelGrid, ...
        pore.connectivity.y.componentIds, 'Pore connected in y', xRange, yRange);

    subplot(2, 2, 3);
    drawConnectivityPanel(xCenters, yCenters, validMask, matrix.mask, matrix.labelGrid, ...
        matrix.connectivity.x.componentIds, 'Matrix connected in x', xRange, yRange);

    subplot(2, 2, 4);
    drawConnectivityPanel(xCenters, yCenters, validMask, matrix.mask, matrix.labelGrid, ...
        matrix.connectivity.y.componentIds, 'Matrix connected in y', xRange, yRange);

    annotation(fig, 'textbox', [0.28 0.95 0.44 0.04], 'String', ...
        sprintf('Directional connectivity highlights @ timestep %g', timestep), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
end

function drawConnectivityPanel(xCenters, yCenters, validMask, phaseMask, labelGrid, componentIds, ttl, xRange, yRange)
    ax = gca;
    highlightMask = false(size(labelGrid));
    if ~isempty(componentIds)
        highlightMask = ismember(labelGrid, componentIds);
    end

    z = nan(size(labelGrid));
    z(validMask & ~phaseMask) = 0;
    z(validMask & phaseMask & ~highlightMask) = 1;
    z(validMask & highlightMask) = 2;

    imagesc(ax, xCenters, yCenters, z);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    applyPlotRanges(ax, xCenters, yCenters, xRange, yRange);
    box(ax, 'on');
    grid(ax, 'off');
    set(ax, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 11);
    colormap(ax, [1.00 1.00 1.00; 0.77 0.84 0.93; 0.88 0.42 0.05]);
    setAxesCLim(ax, [0 2]);
    xlabel(ax, 'x', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax, 'y', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax, ttl, 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    if ~any(highlightMask(:))
        xLimNow = get(ax, 'XLim');
        yLimNow = get(ax, 'YLim');
        text(ax, mean(xLimNow), mean(yLimNow), 'No connected component', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', [0.25 0.25 0.25], 'FontName', 'Times New Roman', ...
            'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1]);
    end
end

function applyPlotRanges(ax, xCenters, yCenters, xRange, yRange)
    if ~isempty(xRange)
        xlim(ax, normalizePlotRange(xRange, xCenters));
    end
    if ~isempty(yRange)
        ylim(ax, normalizePlotRange(yRange, yCenters));
    end
end

function rangeOut = normalizePlotRange(rangeIn, centers)
    rangeOut = double(rangeIn(:).');
    if isempty(rangeOut) || (rangeOut(2) > rangeOut(1))
        return;
    end
    padding = inferPlotPadding(centers);
    rangeOut = [rangeOut(1) - padding, rangeOut(2) + padding];
end

function padding = inferPlotPadding(centers)
    centers = double(centers(:));
    centers = centers(isfinite(centers));
    if numel(centers) >= 2
        d = diff(sort(centers));
        d = d(d > 0);
        if ~isempty(d)
            padding = 0.5 * min(d);
            return;
        end
    end
    padding = 0.5;
end

function fig = plotCountDistribution(dist, phaseName, xRange, histScale, plotStyle)
    fig = [];
    if isempty(dist.count)
        return;
    end
    fig = figure('Color', 'w', 'Name', [phaseName, ' Count Distribution']);
    ax = axes('Parent', fig);
    h = [];
    labels = {};
    hData = plotDistributionSeries(ax, dist.centers, dist.count, ...
        plotStyle, [0.35 0.6 0.85]);
    hold(ax, 'on');
    h(end + 1) = hData;
    labels{end + 1} = 'Count';
    [hFit, fitLabels] = plotFitCurves(ax, dist.fit, 'count');
    h = [h, hFit];
    labels = [labels, fitLabels];
    xlabel(ax, 'Equivalent Diameter', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'Count', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, [phaseName, ' diameter count distribution'], ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    if numel(h) > 1
        legend(ax, h, labels, 'Location', 'best');
    end
    applyHistScale(ax, histScale);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.16, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
    applyDistributionXRange(ax, xRange, dist.centers);
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

function applyHistScale(ax, histScale)
    mode = lower(strtrim(histScale));
    switch mode
        case 'linear'
            set(ax, 'XScale', 'linear', 'YScale', 'linear');
        case 'semilogx'
            set(ax, 'XScale', 'log', 'YScale', 'linear');
        case {'semilogy', 'semilog'}
            set(ax, 'XScale', 'linear', 'YScale', 'log');
        case {'loglog', 'log'}
            set(ax, 'XScale', 'log', 'YScale', 'log');
    end
end

function applyDistributionXRange(ax, xRange, xFallback)
    if isempty(xRange)
        return;
    end
    xlim(ax, normalizePlotRange(xRange, xFallback));
end

function fig = plotPositionDistribution(dist, phaseName, axisName)
    fig = [];
    if isempty(dist.meanEquivDiameter)
        return;
    end
    valid = (dist.count > 0) & isfinite(dist.meanEquivDiameter);
    if ~any(valid)
        return;
    end
    fig = figure('Color', 'w', 'Name', [phaseName, ' Position Distribution ', upper(axisName)]);
    ax = axes('Parent', fig);
    plot(ax, dist.centers(valid), dist.meanEquivDiameter(valid), 'o-', ...
        'LineWidth', 1.5, 'Color', [0.12 0.40 0.75], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.12 0.40 0.75]);
    xlabel(ax, axisName, 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel(ax, 'Mean Equivalent Diameter', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(ax, sprintf('%s mean size vs %s bin', phaseName, axisName), ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.16, 'LineWidth', 1.0, 'FontName', 'Times New Roman', 'FontSize', 12);
end

function fig = plotDirectionalProfile(profile, timestep)
    fig = figure('Color', 'w', 'Name', ['Directional Profile ', upper(profile.axis)]);

    ax1 = subplot(3, 1, 1, 'Parent', fig);
    plot(ax1, profile.centers, profile.porosity, 'o-', ...
        'LineWidth', 1.5, 'Color', [0.10 0.42 0.78], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.10 0.42 0.78]);
    ylabel(ax1, 'Porosity', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax1, sprintf('Pore-structure profile along %s @ timestep %g', profile.axis, timestep), ...
        'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    pd_plot_style_profile_axis(ax1);

    ax2 = subplot(3, 1, 2, 'Parent', fig);
    plot(ax2, profile.centers, profile.specificInterfaceBulk, 's-', ...
        'LineWidth', 1.5, 'Color', [0.85 0.38 0.08], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.85 0.38 0.08]);
    ylabel(ax2, 'Specific Interface', 'FontName', 'Times New Roman', 'FontSize', 12);
    pd_plot_style_profile_axis(ax2);

    ax3 = subplot(3, 1, 3, 'Parent', fig);
    yyaxis(ax3, 'left');
    plot(ax3, profile.centers, profile.connectivityFraction, '^-', ...
        'LineWidth', 1.5, 'Color', [0.15 0.55 0.20], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.15 0.55 0.20]);
    ylabel(ax3, 'Conn. Fraction', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylim(ax3, [0 1]);
    yyaxis(ax3, 'right');
    stairs(ax3, profile.centers, double(profile.connectivityFlag), '--', ...
        'LineWidth', 1.4, 'Color', [0.55 0.10 0.10]);
    ylabel(ax3, 'Connected (0/1)', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylim(ax3, [0 1]);
    xlabel(ax3, profile.axis, 'FontName', 'Times New Roman', 'FontSize', 13);
    legend(ax3, {'Connected area fraction', 'Perp. connected flag'}, 'Location', 'best');
    pd_plot_style_profile_axis(ax3);
end

function plots = emptyPlotsStruct()
    plots = struct( ...
        'phaseFig', [], ...
        'poreLabelFig', [], ...
        'matrixLabelFig', [], ...
        'connectivityFig', [], ...
        'evolutionFig', [], ...
        'profileFigX', [], ...
        'profileFigY', [], ...
        'poreCountFig', [], ...
        'matrixCountFig', [], ...
        'porePositionFigX', [], ...
        'porePositionFigY', [], ...
        'matrixPositionFigX', [], ...
        'matrixPositionFigY', []);
end

function setAxesCLim(ax, limits)
    if nargin < 2 || isempty(limits)
        return;
    end
    set(ax, 'CLim', limits);
end
