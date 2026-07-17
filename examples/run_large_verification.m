function report = run_large_verification(outputDir)
%RUN_LARGE_VERIFICATION Run generated 1D, 2D, cluster, network, and mass-x cases.
% Saves ordinary MATLAB figure dashboards that work on R2016b.

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    postdata_startup();
    addpath(fullfile(rootDir, 'examples'));
    fixtureDir = fullfile(rootDir, 'fixtures', 'generated');
    paths = generate_large_test_data(fixtureDir);
    if nargin < 1 || isempty(outputDir)
        outputDir = fullfile(rootDir, 'outputs', 'verification');
    end
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    oldVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanupObj = onCleanup(@() set(0, 'DefaultFigureVisible', oldVisible)); %#ok<NASGU>
    plotOptions = defaultPlotOptions();

    chunk1d = runRequest('chunk', paths.chunk1d, {}, 'full');
    file1d = fullfile(outputDir, 'large_chunk1d_dashboard.png');
    renderDashboard(chunk1d, plotOptions, file1d);

    massxOptions = {'ChunkDim', '1d', 'InitialDensity', 7.3, ...
        'ParticleSpacing', 0.005, 'RawLengthUnitUm', 10, ...
        'TransverseWidth', 2.46, 'CoordinateFactor', 10};
    massx = runRequest('massx', paths.chunk1d, massxOptions, 'full');
    if ~all(massx.isMonotonicDecreasing)
        error('postdata:LargeVerificationMassXMonotonicity', ...
            'Large mass-x cumulative curves must be monotonic decreasing.');
    end
    fileMassX = fullfile(outputDir, 'large_massx_dashboard.png');
    renderDashboard(massx, plotOptions, fileMassX);

    massx2dOptions = {'ChunkDim', '2d', 'InitialDensity', 7.3, ...
        'ParticleSpacing', 0.01, 'RawLengthUnitUm', 10, ...
        'CoordinateFactor', 10, 'SliceCentersY', [0.6 1.2 1.8], ...
        'SliceWidthsY', 0.6};
    massx2d = runRequest('massx', paths.chunk2d, massx2dOptions, 'full');
    if ~all(massx2d.isMonotonicDecreasing)
        error('postdata:LargeVerificationMassX2dMonotonicity', ...
            'Large sliced 2D mass-x curves must be monotonic decreasing.');
    end
    fileMassX2d = fullfile(outputDir, 'large_massx2d_dashboard.png');
    renderDashboard(massx2d, plotOptions, fileMassX2d);

    massv = runRequest('mass-v', paths.massv, {}, 'full');
    if ~all(massv.isMonotonicDecreasing)
        error('postdata:LargeVerificationMassVMonotonicity', ...
            'Default large mass-v cumulative curves must be monotonic decreasing.');
    end
    fileMassV = fullfile(outputDir, 'large_massv_dashboard.png');
    renderDashboard(massv, plotOptions, fileMassV);

    chunk2d = runRequest('chunk', paths.chunk2d, {}, 'full');
    file2d = fullfile(outputDir, 'large_chunk2d_dashboard.png');
    renderDashboard(chunk2d, plotOptions, file2d);

    clusterOptions = {'Dim', 2, 'Dx', 0.025, 'MeanNumBins', 30, ...
        'DiameterHistBinSize', 0.05, 'XVarForMean', 'c_x'};
    cluster = runRequest('cluster', paths.cluster, clusterOptions, 'full');
    fileCluster = fullfile(outputDir, 'large_cluster_dashboard.png');
    renderDashboard(cluster, plotOptions, fileCluster);

    networkOptions = {'ThresholdN', 1, 'GeometryMode', 'original', ...
        'PositionAxis', 'both', 'ProfileAxis', 'both', ...
        'PositionNumBins', 30, 'ProfileNumBins', 30};
    network = runRequest('network2d', paths.chunk2d, networkOptions, 'full');
    fileNetwork = fullfile(outputDir, 'large_network2d_dashboard.png');
    renderDashboard(network, plotOptions, fileNetwork);

    report = struct();
    report.files = struct('chunk1d', file1d, 'massx', fileMassX, ...
        'massx2d', fileMassX2d, 'massv', fileMassV, ...
        'chunk2d', file2d, 'cluster', fileCluster, 'network2d', fileNetwork);
    finite1d = chunk1d.value(isfinite(chunk1d.value));
    report.chunk1d = struct('rows', numel(chunk1d.x), ...
        'minimum', min(finite1d), 'maximum', max(finite1d));
    report.massx = struct('rows', numel(massx.coordinate), ...
        'totalArealDensity', massx.cumulativeDensity(1, end), ...
        'isMonotonicDecreasing', all(massx.isMonotonicDecreasing));
    report.massx2d = struct('rows', numel(massx2d.coordinate), ...
        'slices', size(massx2d.cumulativeDensity, 2), ...
        'totalArealDensity', massx2d.cumulativeDensity(1, :), ...
        'isMonotonicDecreasing', all(massx2d.isMonotonicDecreasing));
    report.massv = struct('rows', numel(massv.velocity), ...
        'totalArealDensity', massv.cumulativeDensity(1, end), ...
        'isMonotonicDecreasing', all(massv.isMonotonicDecreasing));
    finite2d = chunk2d.value(isfinite(chunk2d.value));
    report.chunk2d = struct('rows', numel(chunk2d.x), ...
        'minimum', min(finite2d), 'maximum', max(finite2d));
    report.cluster = struct('selectedRows', cluster.selectedRows, ...
        'meanDiameter', mean(cluster.diameter), 'maximumDiameter', max(cluster.diameter));
    report.network2d = struct('validCells', network.global.numValidCells, ...
        'porosity', network.global.porosity, ...
        'poreComponents', network.pore.numComponents, ...
        'matrixComponents', network.matrix.numComponents);
    writeReport(fullfile(outputDir, 'large_verification_summary.txt'), report);
end

function result = runRequest(type, filePath, analysisOptions, resultLevel)
    request = pd_create_request(type);
    request.filePath = filePath;
    request.makePlots = false;
    request.progressMode = 'off';
    request.resultLevel = resultLevel;
    request.analysisOptions = analysisOptions;
    result = postdata_run(request);
end

function options = defaultPlotOptions()
    catalog = pd_plot_option_catalog();
    options = pd_plot_options_from_table(catalog, pd_catalog_table_data(catalog));
    options.UpdateMode = 'replace';
end

function renderDashboard(result, options, filePath)
    views = pd_result_plot_views(result);
    count = numel(views);
    [rows, columns] = pd_subplot_grid(count);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [50, 50, 520 * columns, 390 * rows]);
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    gapX = 0.07;
    gapY = 0.12;
    marginX = 0.07;
    marginY = 0.12;
    width = (1 - 2 * marginX - (columns - 1) * gapX) / columns;
    height = (1 - 2 * marginY - (rows - 1) * gapY) / rows;
    for i = 1:count
        row = floor((i - 1) / columns);
        column = mod(i - 1, columns);
        left = marginX + column * (width + gapX);
        bottom = 1 - marginY - (row + 1) * height - row * gapY;
        ax = axes('Parent', fig, 'Units', 'normalized', ...
            'Position', [left, bottom, width, height]);
        pd_render_result(ax, result, options, views(i).Id);
    end
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, filePath, '-dpng', '-r130');
end

function writeReport(filePath, report)
    fid = fopen(filePath, 'w');
    if fid < 0, error('postdata:ReportOpenFailed', 'Cannot create %s.', filePath); end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'POST_DATA2 large verification\n');
    fprintf(fid, 'chunk1d.rows=%d\n', report.chunk1d.rows);
    fprintf(fid, 'chunk1d.range=[%.8g, %.8g]\n', ...
        report.chunk1d.minimum, report.chunk1d.maximum);
    fprintf(fid, 'massx.rows=%d\n', report.massx.rows);
    fprintf(fid, 'massx.totalArealDensity=%.8g mg/cm^2\n', ...
        report.massx.totalArealDensity);
    fprintf(fid, 'massx.isMonotonicDecreasing=%d\n', ...
        report.massx.isMonotonicDecreasing);
    fprintf(fid, 'massx2d.rows=%d\n', report.massx2d.rows);
    fprintf(fid, 'massx2d.slices=%d\n', report.massx2d.slices);
    fprintf(fid, 'massx2d.totalArealDensity=%s mg/cm^2\n', ...
        mat2str(report.massx2d.totalArealDensity, 8));
    fprintf(fid, 'massx2d.isMonotonicDecreasing=%d\n', ...
        report.massx2d.isMonotonicDecreasing);
    fprintf(fid, 'massv.rows=%d\n', report.massv.rows);
    fprintf(fid, 'massv.totalArealDensity=%.8g mg/cm^2\n', ...
        report.massv.totalArealDensity);
    fprintf(fid, 'massv.isMonotonicDecreasing=%d\n', ...
        report.massv.isMonotonicDecreasing);
    fprintf(fid, 'chunk2d.rows=%d\n', report.chunk2d.rows);
    fprintf(fid, 'chunk2d.range=[%.8g, %.8g]\n', ...
        report.chunk2d.minimum, report.chunk2d.maximum);
    fprintf(fid, 'cluster.selectedRows=%d\n', report.cluster.selectedRows);
    fprintf(fid, 'cluster.meanDiameter=%.8g\n', report.cluster.meanDiameter);
    fprintf(fid, 'cluster.maximumDiameter=%.8g\n', report.cluster.maximumDiameter);
    fprintf(fid, 'network2d.validCells=%d\n', report.network2d.validCells);
    fprintf(fid, 'network2d.porosity=%.8g\n', report.network2d.porosity);
    fprintf(fid, 'network2d.poreComponents=%d\n', report.network2d.poreComponents);
    fprintf(fid, 'network2d.matrixComponents=%d\n', report.network2d.matrixComponents);
end
