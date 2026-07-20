function files = pd_export_result(result, ax, options)
%PD_EXPORT_RESULT Export result data and the current visualization.
% Supports MAT, summary/detail CSV, PNG, FIG, PDF, and manifest output.

    if nargin < 2, ax = []; end
    if nargin < 3, options = struct(); end
    options = normalizeOutputOptions(options);
    outputDir = pd_to_char(options.Directory);
    if ~exist(outputDir, 'dir')
        [ok, message] = mkdir(outputDir);
        if ~ok
            error('postdata:CreateOutputDirectoryFailed', '%s', message);
        end
    end

    prefix = regexprep(pd_to_char(options.Prefix), '[^A-Za-z0-9_-]', '_');
    if isempty(prefix)
        prefix = 'postdata';
    end
    stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<DATST,TNOW1>
    baseName = sprintf('%s_%s_%s', prefix, result.analysisType, stamp);
    basePath = uniqueBasePath(fullfile(outputDir, baseName));
    [~, committedBaseName] = fileparts(basePath);
    stageDir = tempname(outputDir);
    [ok, message] = mkdir(stageDir);
    if ~ok
        error('postdata:CreateExportStageFailed', '%s', message);
    end
    stageCleanup = onCleanup(@() removeStageDirectory(stageDir)); %#ok<NASGU>
    stageBasePath = fullfile(stageDir, committedBaseName);
    stagedFiles = {};

    if options.SaveMAT
        matPath = [stageBasePath, '.mat'];
        save(matPath, 'result');
        stagedFiles{end + 1} = matPath; %#ok<AGROW>
    end
    if options.SaveCSV
        csvPath = [stageBasePath, '.csv'];
        writeSummaryCsv(csvPath, pd_result_summary(result));
        stagedFiles{end + 1} = csvPath; %#ok<AGROW>
    end
    if options.SaveDetailCSV
        detailFiles = pd_export_detail_csv(result, stageBasePath);
        stagedFiles = [stagedFiles, detailFiles]; %#ok<AGROW>
    end
    if options.SavePlotDataCSV
        plotDataFiles = pd_export_plot_data_csv(result, stageBasePath);
        stagedFiles = [stagedFiles, plotDataFiles]; %#ok<AGROW>
    end

    needsFigure = options.SavePNG || options.SaveFIG || options.SavePDF;
    exportFigure = [];
    if needsFigure
        if nargin < 2 || isempty(ax) || ~all(ishghandle(ax))
            error('postdata:MissingExportAxes', 'Valid axes are required for figure export.');
        end
        exportFigure = makeExportFigure(ax, options);
        cleanupFigure = onCleanup(@() close(exportFigure)); %#ok<NASGU>
    end
    if options.SavePNG
        pngPath = [stageBasePath, '.png'];
        print(exportFigure, pngPath, '-dpng', ['-r', num2str(round(options.DPI))]);
        stagedFiles{end + 1} = pngPath; %#ok<AGROW>
    end
    if options.SaveFIG
        figPath = [stageBasePath, '.fig'];
        savefig(exportFigure, figPath);
        stagedFiles{end + 1} = figPath; %#ok<AGROW>
    end
    if options.SavePDF
        pdfPath = [stageBasePath, '.pdf'];
        print(exportFigure, pdfPath, '-dpdf');
        stagedFiles{end + 1} = pdfPath; %#ok<AGROW>
    end
    if options.SaveManifest
        manifestPath = [stageBasePath, '_manifest.txt'];
        pd_write_manifest(manifestPath, result, options);
        stagedFiles{end + 1} = manifestPath; %#ok<AGROW>
    end
    files = commitStagedFiles(stagedFiles, outputDir);
end

function options = normalizeOutputOptions(input)
    if ~isstruct(input) || ~isscalar(input)
        error('postdata:BadOutputOptions', ...
            'Output options must be a scalar structure.');
    end
    catalog = pd_output_option_catalog();
    data = pd_catalog_table_data(catalog);
    names = fieldnames(input);
    for i = 1:numel(names)
        row = find(strcmpi(names{i}, {catalog.Name}), 1, 'first');
        if isempty(row)
            error('postdata:UnknownOutputOption', ...
                'Unknown output option: %s.', names{i});
        end
        data{row, 2} = input.(names{i});
    end
    options = pd_output_options_from_table(catalog, data);
end

function files = commitStagedFiles(stagedFiles, outputDir)
    files = cell(size(stagedFiles));
    committed = {};
    try
        for i = 1:numel(stagedFiles)
            [~, name, extension] = fileparts(stagedFiles{i});
            destination = fullfile(outputDir, [name, extension]);
            [ok, message] = movefile(stagedFiles{i}, destination);
            if ~ok
                error('postdata:CommitExportFailed', '%s', message);
            end
            files{i} = destination;
            committed{end + 1} = destination; %#ok<AGROW>
        end
    catch err
        for i = 1:numel(committed)
            if exist(committed{i}, 'file')
                try
                    delete(committed{i});
                catch
                    % Preserve the original export error.
                end
            end
        end
        rethrow(err);
    end
end

function removeStageDirectory(stageDir)
    if exist(stageDir, 'dir')
        try
            rmdir(stageDir, 's');
        catch
            % Cleanup failure must not hide the export result or root cause.
        end
    end
end

function writeSummaryCsv(filePath, data)
    fid = fopen(filePath, 'w');
    if fid < 0
        error('postdata:CsvOpenFailed', 'Cannot create CSV: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'Field,Value\n');
    for i = 1:size(data, 1)
        fprintf(fid, '%s,%s\n', csvEscape(data{i, 1}), csvEscape(data{i, 2}));
    end
end

function text = csvEscape(value)
    text = pd_to_char(value);
    text = strrep(text, '"', '""');
    text = ['"', text, '"'];
end

function fig = makeExportFigure(ax, options)
    widthCm = options.FigureWidthCm;
    heightCm = options.FigureHeightCm;
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, widthCm, heightCm]);
    copiedAxes = copyobj(ax, fig);
    count = numel(copiedAxes);
    [rows, columns] = pd_subplot_grid(count);
    gapX = 0.07;
    gapY = 0.09;
    marginX = 0.08;
    marginY = 0.09;
    width = (1 - 2 * marginX - (columns - 1) * gapX) / columns;
    height = (1 - 2 * marginY - (rows - 1) * gapY) / rows;
    for i = 1:count
        row = floor((i - 1) / columns);
        column = mod(i - 1, columns);
        left = marginX + column * (width + gapX);
        bottom = 1 - marginY - (row + 1) * height - row * gapY;
        set(copiedAxes(i), 'Units', 'normalized', ...
            'Position', [left, bottom, width, height]);
    end
    set(fig, 'PaperUnits', 'centimeters', ...
        'PaperPositionMode', 'manual', ...
        'PaperPosition', [0, 0, widthCm, heightCm], ...
        'PaperSize', [widthCm, heightCm]);
end

function output = uniqueBasePath(basePath)
    output = basePath;
    index = 1;
    while ~isempty(dir([output, '*']))
        output = sprintf('%s_%03d', basePath, index);
        index = index + 1;
    end
end
