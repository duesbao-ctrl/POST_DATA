function out = demo_network2d_periodic_topology()
% demo_network2d_periodic_topology
% -------------------------------------------------------------
% Build a deterministic periodic 2D pore network, run the paper-grade
% network2d analysis, and save the results for inspection.

    repoDir = fileparts(mfilename('fullpath'));
    outputDir = fullfile(repoDir, 'tmp_network2d_periodic_topology_demo');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    rng(7);
    nx = 28;
    ny = 24;
    generationFieldThreshold = 0.25;
    analysisThresholdN = 1;

    [xGrid, yGrid] = meshgrid((0:nx-1) / nx, (0:ny-1) / ny);
    field = buildPeriodicField(xGrid, yGrid);
    poreMask = field > generationFieldThreshold;

    ncountGrid = 2 * ones(ny, nx);
    ncountGrid(poreMask) = 0;

    chunkFile = fullfile(outputDir, 'bin2d_periodic_complex_demo_dx_1_dy_1.txt');
    writeChunk2dCase(chunkFile, ncountGrid, 100);

    out = run_analysis('network2d', ...
        'BaseDir', outputDir, ...
        'ChunkDim', '2d', ...
        'ChunkFile', 'bin2d_periodic_complex_demo_dx_1_dy_1.txt', ...
        'SelectBy', 'Index', ...
        'Index', 1, ...
        'ProgressMode', 'off', ...
        'NetworkOptions', {'ThresholdN', analysisThresholdN, ...
                           'Boundary', 'periodic-xy', ...
                           'EnableSkeletonGraph', true, ...
                           'MeanPowerM', 1, 'MeanPowerN', 0, ...
                           'PositionAxis', 'both', ...
                           'ProfileAxis', 'both', ...
                           'PositionNumBins', 12, ...
                           'ProfileNumBins', 12, ...
                           'MakePlots', false});

    save(fullfile(outputDir, 'network2d_periodic_complex_demo_results.mat'), 'out');
    writeTextFile(fullfile(outputDir, 'network2d_periodic_complex_demo_stats.json'), ...
        jsonencode(out.stats, 'PrettyPrint', true));
    writeTextFile(fullfile(outputDir, 'network2d_periodic_complex_demo_stats.txt'), ...
        evalc('disp(out.stats)'));
    writeTextFile(fullfile(outputDir, 'network2d_periodic_complex_demo_summary.txt'), ...
        buildSummaryText(out, chunkFile, outputDir, generationFieldThreshold, analysisThresholdN));
    savePhasePreview(fullfile(outputDir, 'network2d_periodic_complex_demo_phase.png'), poreMask);

    printSummaryLines(buildConsoleSummary(out, chunkFile, outputDir, generationFieldThreshold, analysisThresholdN));
end

function field = buildPeriodicField(xGrid, yGrid)
    field = zeros(size(xGrid));
    for ax = 0:4
        for ay = 0:4
            if ax == 0 && ay == 0
                continue;
            end
            amp = randn() / (1 + ax^2 + ay^2)^0.7;
            phase1 = 2 * pi * rand();
            phase2 = 2 * pi * rand();
            field = field + amp * cos(2 * pi * (ax * xGrid + ay * yGrid) + phase1);
            field = field + 0.6 * amp * sin(2 * pi * (ax * xGrid - ay * yGrid) + phase2);
        end
    end

    field = field ./ std(field(:));
end

function writeChunk2dCase(filePath, ncountGrid, timestep)
    fid = fopen(filePath, 'w');
    if fid < 0
        error('demo_network2d_periodic_topology:FileOpenFail', ...
            'Cannot open %s for writing.', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '# Chunk-averaged data for fix temp/chunk and group all\n');
    fprintf(fid, '# Timestep Number-of-chunks Total-count\n');
    fprintf(fid, '# Chunk Coord1 Coord2 Ncount\n');
    fprintf(fid, '%g %d %g\n', timestep, numel(ncountGrid), sum(ncountGrid(:)));

    chunkId = 1;
    for iy = 1:size(ncountGrid, 1)
        for ix = 1:size(ncountGrid, 2)
            fprintf(fid, '%d %g %g %g\n', chunkId, ix - 1, iy - 1, ncountGrid(iy, ix));
            chunkId = chunkId + 1;
        end
    end
end

function savePhasePreview(filePath, poreMask)
    fig = figure('Visible', 'off', 'Color', 'w');
    imagesc(poreMask);
    axis image tight;
    set(gca, 'YDir', 'normal');
    colormap(gca, [0.11 0.34 0.62; 0.88 0.34 0.03]);
    title('Periodic-XY Complex Pore Network');
    xlabel('x cell index');
    ylabel('y cell index');
    exportgraphics(fig, filePath, 'Resolution', 200);
    close(fig);
end

function writeTextFile(filePath, textContent)
    fid = fopen(filePath, 'w');
    if fid < 0
        error('demo_network2d_periodic_topology:FileOpenFail', ...
            'Cannot open %s for writing.', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid, textContent, 'char');
end

function textOut = buildSummaryText(out, chunkFile, outputDir, generationFieldThreshold, analysisThresholdN)
    lines = buildConsoleSummary(out, chunkFile, outputDir, generationFieldThreshold, analysisThresholdN);
    textOut = strjoin(lines, newline);
end

function lines = buildConsoleSummary(out, chunkFile, outputDir, generationFieldThreshold, analysisThresholdN)
    lines = {
        'network2d periodic topology demo'
        ['chunkFile: ', chunkFile]
        ['outputDir: ', outputDir]
        ['generationFieldThreshold: ', num2str(generationFieldThreshold)]
        ['analysisThresholdN: ', num2str(analysisThresholdN)]
        ['boundary: ', out.stats.meta.boundary]
        ['phi: ', num2str(out.stats.geometry.phi, '%.6f')]
        ['interfaceLength: ', num2str(out.stats.geometry.interfaceLength, '%.6f')]
        ['specificInterface: ', num2str(out.stats.geometry.specificInterface, '%.6f')]
        ['pore beta0/beta1/chi: ', sprintf('%g / %g / %g', ...
            out.stats.topology.pore.beta0, out.stats.topology.pore.beta1, out.stats.topology.pore.chi)]
        ['pore largestFraction: ', num2str(out.stats.connectivity.pore.largestFraction, '%.6f')]
        ['pore percolatesX/percolatesY: ', sprintf('%d / %d', ...
            logical(out.stats.connectivity.pore.percolatesX), logical(out.stats.connectivity.pore.percolatesY))]
        ['pore wrapsX/wrapsY: ', sprintf('%d / %d', logical(out.pore.wrapsX), logical(out.pore.wrapsY))]
        ['pore branchPoints/endPoints: ', sprintf('%g / %g', ...
            out.stats.network.pore.branchPoints, out.stats.network.pore.endPoints)]
        ['matrix beta0/beta1/chi: ', sprintf('%g / %g / %g', ...
            out.stats.topology.matrix.beta0, out.stats.topology.matrix.beta1, out.stats.topology.matrix.chi)]
        ['matrix largestFraction: ', num2str(out.stats.connectivity.matrix.largestFraction, '%.6f')]
        ['matrix percolatesX/percolatesY: ', sprintf('%d / %d', ...
            logical(out.stats.connectivity.matrix.percolatesX), logical(out.stats.connectivity.matrix.percolatesY))]
        ['matrix wrapsX/wrapsY: ', sprintf('%d / %d', logical(out.matrix.wrapsX), logical(out.matrix.wrapsY))]
        ['matrix thickness min/p5/mean: ', sprintf('%.6f / %.6f / %.6f', ...
            out.stats.thickness.matrix.min, out.stats.thickness.matrix.p5, out.stats.thickness.matrix.mean)]
        ['matrix fragmentation count/largestFraction: ', sprintf('%g / %.6f', ...
            out.stats.fragmentation.matrix.count, out.stats.fragmentation.matrix.largestFraction)]
        ['statsJson: ', fullfile(outputDir, 'network2d_periodic_complex_demo_stats.json')]
        ['statsTxt: ', fullfile(outputDir, 'network2d_periodic_complex_demo_stats.txt')]
        ['resultsMat: ', fullfile(outputDir, 'network2d_periodic_complex_demo_results.mat')]
        ['phasePng: ', fullfile(outputDir, 'network2d_periodic_complex_demo_phase.png')]
        };
end

function printSummaryLines(lines)
    for i = 1:numel(lines)
        fprintf('%s\n', lines{i});
    end
end
