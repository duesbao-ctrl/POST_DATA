function colors = pd_resolve_plot_palette(options)
%PD_RESOLVE_PLOT_PALETTE Resolve a named or custom RGB series palette.

    if isfield(options, 'CustomColorOrder') && ~isempty(options.CustomColorOrder)
        colors = options.CustomColorOrder;
    else
        colors = pd_plot_palette(options.ColorPalette);
    end
end
