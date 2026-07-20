function fig = pd_render_network2d_evolution(evolution)
    fig = [];
    if isempty(evolution.axisValue)
        return;
    end

    x = evolution.axisValue;
    if any(~isfinite(x))
        x = evolution.stepIndex;
        xLabel = 'index';
    else
        xLabel = evolution.axis;
    end

    fig = figure('Color', 'w', 'Name', 'Network2d Evolution');

    ax1 = subplot(2, 2, 1, 'Parent', fig);
    plot(ax1, x, evolution.geometry.phi, 'o-', 'LineWidth', 1.5, ...
        'Color', [0.10 0.42 0.78], 'MarkerFaceColor', [0.10 0.42 0.78]);
    xlabel(ax1, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax1, '\phi', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax1, 'Porosity Evolution', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    pd_plot_style_profile_axis(ax1);

    ax2 = subplot(2, 2, 2, 'Parent', fig);
    plot(ax2, x, evolution.geometry.specificInterface, 's-', 'LineWidth', 1.5, ...
        'Color', [0.85 0.38 0.08], 'MarkerFaceColor', [0.85 0.38 0.08]);
    xlabel(ax2, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax2, 'L / area', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax2, 'Specific Interface Evolution', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    pd_plot_style_profile_axis(ax2);

    ax3 = subplot(2, 2, 3, 'Parent', fig);
    plot(ax3, x, evolution.connectivity.pore.largestFraction, 'o-', 'LineWidth', 1.5, ...
        'Color', [0.15 0.55 0.20], 'MarkerFaceColor', [0.15 0.55 0.20]);
    hold(ax3, 'on');
    plot(ax3, x, evolution.connectivity.matrix.largestFraction, 'd-', 'LineWidth', 1.5, ...
        'Color', [0.55 0.10 0.10], 'MarkerFaceColor', [0.55 0.10 0.10]);
    xlabel(ax3, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax3, 'Largest Fraction', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax3, 'Largest Cluster Fraction', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    legend(ax3, {'Pore', 'Matrix'}, 'Location', 'best');
    pd_plot_style_profile_axis(ax3);

    ax4 = subplot(2, 2, 4, 'Parent', fig);
    plot(ax4, x, evolution.topology.pore.beta0, 'o-', 'LineWidth', 1.3, ...
        'Color', [0.10 0.42 0.78], 'MarkerFaceColor', [0.10 0.42 0.78]);
    hold(ax4, 'on');
    plot(ax4, x, evolution.topology.pore.beta1, 's-', 'LineWidth', 1.3, ...
        'Color', [0.85 0.38 0.08], 'MarkerFaceColor', [0.85 0.38 0.08]);
    plot(ax4, x, evolution.topology.matrix.beta0, 'd--', 'LineWidth', 1.3, ...
        'Color', [0.15 0.55 0.20], 'MarkerFaceColor', [0.15 0.55 0.20]);
    plot(ax4, x, evolution.topology.matrix.beta1, '^-', 'LineWidth', 1.3, ...
        'Color', [0.55 0.10 0.10], 'MarkerFaceColor', [0.55 0.10 0.10]);
    xlabel(ax4, xLabel, 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel(ax4, 'Betti Number', 'FontName', 'Times New Roman', 'FontSize', 12);
    title(ax4, 'Topology Evolution', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    legend(ax4, {'pore \beta_0', 'pore \beta_1', 'matrix \beta_0', 'matrix \beta_1'}, 'Location', 'best');
    pd_plot_style_profile_axis(ax4);
end
