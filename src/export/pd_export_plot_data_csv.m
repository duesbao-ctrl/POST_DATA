function files = pd_export_plot_data_csv(result, basePath)
%PD_EXPORT_PLOT_DATA_CSV Export the exact numeric table behind every view.

    views = pd_result_plot_views(result);
    files = cell(1, numel(views));
    for i = 1:numel(views)
        safeId = regexprep(views(i).Id, '[^A-Za-z0-9_-]', '_');
        files{i} = [basePath, '_view_', safeId, '.csv'];
        pd_write_plot_data_csv(files{i}, ...
            pd_result_plot_data(result, views(i).Id));
    end
end
