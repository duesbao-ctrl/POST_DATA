function pd_write_plot_data_csv(filePath, plotData, rowIndices)
%PD_WRITE_PLOT_DATA_CSV Write one plot-view data contract to CSV.

    if nargin < 3
        rowIndices = [];
    end
    text = pd_plot_data_text(plotData, ',', rowIndices);
    fid = fopen(filePath, 'w');
    if fid < 0
        error('postdata:CsvOpenFailed', 'Cannot create CSV: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', text);
end
