function colors = pd_plot_palette(name)
%PD_PLOT_PALETTE Return a publication-friendly, R2016b-compatible palette.

    switch lower(strtrim(pd_to_char(name)))
        case 'colorblind'
            colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
                0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; ...
                0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; ...
                0.6350 0.0780 0.1840];
        case 'grayscale'
            colors = [0.10 0.10 0.10; 0.35 0.35 0.35; ...
                0.58 0.58 0.58; 0.78 0.78 0.78];
        case 'highcontrast'
            colors = [0 0 0; 0.80 0 0; 0 0.35 0.75; ...
                0 0.55 0.20; 0.75 0.45 0];
        otherwise
            colors = lines(7);
    end
end
