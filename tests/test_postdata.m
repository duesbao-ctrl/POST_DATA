function tests = test_postdata
%TEST_POSTDATA Integration tests for the single POST_DATA implementation.
% MATLAB R2016b compatible.

    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    postdata_startup();
    testCase.TestData.RootDir = rootDir;
    testCase.TestData.FixtureDir = fullfile(rootDir, 'fixtures', 'sample');
end

function testVersion(testCase)
    info = pd_version();
    verifyEqual(testCase, info.minimumMatlabRelease, 'R2016b');
    verifyEqual(testCase, info.resultSchemaVersion, 1);
end

function testUiTextUsesExplicitUtf8Resource(testCase)
    notRunYet = pd_ui_text('Not run yet');
    verifyEqual(testCase, double(notRunYet), [23578 26410 36816 34892]);
    totalViews = pd_ui_text('Total: %d views', 2);
    verifyEqual(testCase, double(totalViews), [20849 32 50 32 24133]);
    verifyEqual(testCase, pd_ui_text('Unmapped ASCII fallback'), ...
        'Unmapped ASCII fallback');
end

function testExecutableMatlabSourcesAreAsciiOnly(testCase)
    folders = regexp(genpath(testCase.TestData.RootDir), pathsep, 'split');
    checked = 0;
    for folderIndex = 1:numel(folders)
        folder = folders{folderIndex};
        if isempty(folder), continue; end
        files = dir(fullfile(folder, '*.m'));
        for fileIndex = 1:numel(files)
            filePath = fullfile(folder, files(fileIndex).name);
            fid = fopen(filePath, 'rb');
            verifyGreaterThanOrEqual(testCase, fid, 0, filePath);
            cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
            bytes = fread(fid, Inf, '*uint8');
            verifyTrue(testCase, all(bytes <= 127), filePath);
            clear cleanupObj;
            checked = checked + 1;
        end
    end
    verifyGreaterThan(testCase, checked, 0);
end

function testLegacyRunAnalysisUsesModularCore(testCase)
    result = run_analysis('mass-x', 'BaseDir', testCase.TestData.FixtureDir, ...
        'MassXFile', 'bin1d_dx_0.5.txt', ...
        'MassXOptions', [sampleMassXOptions(), {'MakePlot', false}], ...
        'ProgressMode', 'off');
    verifyEqual(testCase, result.analysisType, 'massx');
    verifyEqual(testCase, result.cumulativeDirection, 'high-to-low');

    result = run_analysis('chunk', 'BaseDir', testCase.TestData.FixtureDir, ...
        'ChunkFile', 'bin1d_dx_0.5.txt', 'ChunkDim', 'auto', ...
        'Variable', 'c_rho', 'DoPlot', false, 'ProgressMode', 'off');
    verifyEqual(testCase, result.analysisType, 'chunk');
    verifyEqual(testCase, result.dimension, '1d');
end

function testLegacyWrapperUsesEmbeddedSpidTimeWithoutSlurm(testCase)
    result = run_analysis('chunk', ...
        'BaseDir', testCase.TestData.FixtureDir, ...
        'ChunkFile', 'spid_spatial3d.txt', ...
        'ChunkDim', 'auto', 'Variable', 'rho', ...
        'SelectBy', 'Time', 'Time', 0.19, ...
        'DoPlot', false, 'ProgressMode', 'off');
    verifyEqual(testCase, result.timestep, 200);
    verifyEqual(testCase, result.physicalTime, 0.2, 'AbsTol', 1e-12);
end

function testGridConstruction(testCase)
    [grid, valid, xc, yc, count] = pd_network_build_grid( ...
        [0; 1; 0; 1], [0; 0; 1; 1], [1; 2; 3; 4]);
    verifyEqual(testCase, grid, [1, 2; 3, 4]);
    verifyTrue(testCase, all(valid(:)));
    verifyEqual(testCase, xc, [0, 1]);
    verifyEqual(testCase, yc, [0, 1]);
    verifyEqual(testCase, count, 4);
end

function testDuplicateGridRejected(testCase)
    f = @() pd_network_build_grid([0; 0], [0; 0], [1; 2]);
    verifyError(testCase, f, 'analyze_chunk_network2d:DuplicateGridCell');
end

function testChunkAnalysis(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.chunk.variable = 'c_rho';
    result = postdata_run(request);
    verifyEqual(testCase, result.analysisType, 'chunk');
    verifyEqual(testCase, sum(result.value(isfinite(result.value))), 3.7, 'AbsTol', 1e-12);
end

function testSphRawFieldBypassesMdVolumeDerivation(testCase)
    request = makeRequest(testCase, 'chunk', ...
        'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.chunk.variable = 'c_rho';
    request.chunk.dV = [];
    rawResult = postdata_run(request);
    verifyEqual(testCase, rawResult.sourceType, 'raw');

    request.chunk.variable = 'pressure';
    verifyError(testCase, @() postdata_run(request), ...
        'analyze_chunk_field:NeedDV');
end

function testChunk1dAutoDimensionDoesNotRequireCoord2(testCase)
    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    request.chunk.dimension = 'auto';
    request.chunk.variable = 'c_rho';
    result = postdata_run(request);
    verifyEqual(testCase, result.dimension, '1d');
    verifyEqual(testCase, result.dimensionRequested, 'auto');
    verifyEmpty(testCase, result.y);
    verifyEqual(testCase, result.x, (0:0.5:2).', 'AbsTol', 1e-12);
end

function testCurrentSpidPreambleAndEmbeddedTime(testCase)
    filePath = fullfile(testCase.TestData.FixtureDir, 'spid_spatial3d.txt');
    step = read_chunk_step_fast(filePath, 'SelectBy', 'Time', ...
        'Time', 0.19, 'ProgressMode', 'off');
    verifyEqual(testCase, step.timestep, 200);
    verifyEqual(testCase, step.physicalTime, 0.2, 'AbsTol', 1e-12);
    verifyEqual(testCase, step.inputFormat, 'spid-chunk');
    verifyEqual(testCase, step.unitSystem, 'microscale');
    verifyEqual(testCase, step.taskName, 'bin3d');
    verifyEqual(testCase, step.chunkKind, 'spatial');
    verifyEqual(testCase, numel(step.headerLines), 5);

    full = read_bin_chunk(filePath, 'ProgressMode', 'off');
    verifyEqual(testCase, numel(full.steps), 2);
    verifyEqual(testCase, [full.steps.physicalTime], [0.1 0.2], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, full.metadata.kind, 'spatial');
end

function testCurrentSpidSpatial3dAnalysis(testCase)
    request = makeRequest(testCase, 'chunk', 'spid_spatial3d.txt');
    request.chunk.variable = 'rho';
    result = postdata_run(request);
    verifyEqual(testCase, result.dimension, '3d');
    verifyEqual(testCase, result.coordZ, [0;0;0;0;1;1;1;1]);
    verifyEqual(testCase, result.physicalTime, 0.1, 'AbsTol', 1e-12);
    data = pd_result_plot_data(result, 'field');
    verifyEqual(testCase, data.Kind, 'field3d');
    verifyEqual(testCase, data.ColumnNames, {'x','y','z','rho'});
    verifyEqual(testCase, size(data.Values), [8 4]);
end

function testCurrentSpidPhysicalTimeRequestNeedsNoSlurm(testCase)
    request = makeRequest(testCase, 'chunk', 'spid_spatial3d.txt');
    request.chunk.variable = 'rho';
    request.selection.mode = 'Time';
    request.selection.value = 0.19;
    request.selection.slurmPath = '';
    result = postdata_run(request);
    verifyEqual(testCase, result.timestep, 200);
    verifyEqual(testCase, result.physicalTime, 0.2, 'AbsTol', 1e-12);
end

function testLegacyPhysicalTimeRequestStillRequiresSlurm(testCase)
    folder = tempname;
    mkdir(folder);
    cleanupObj = onCleanup(@() cleanupDirectory(folder)); %#ok<NASGU>
    copyfile(fullfile(testCase.TestData.FixtureDir, ...
        'bin1d_dx_0.5.txt'), folder);
    request = pd_create_request('chunk');
    request.baseDir = folder;
    request.filePath = 'bin1d_dx_0.5.txt';
    request.progressMode = 'off';
    request.makePlots = false;
    request.chunk.variable = 'c_rho';
    request.selection.mode = 'Time';
    request.selection.value = 0.1;
    request.selection.slurmPath = '';
    verifyError(testCase, @() postdata_run(request), ...
        'read_chunk_step_fast:MissingSlurmPath');
end

function testRelativeSlurmPathResolvesAgainstBaseDirectory(testCase)
    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    request.chunk.variable = 'c_rho';
    request.selection.mode = 'Time';
    request.selection.value = 15;
    request.selection.slurmPath = 'slurm_case.log';
    validated = pd_validate_request(request);
    verifyEqual(testCase, validated.selection.slurmPath, ...
        fullfile(testCase.TestData.FixtureDir, 'slurm_case.log'));
    result = postdata_run(request);
    verifyTrue(testCase, any(result.timestep == [100 200]));
end

function testMalformedSpidPhysicalTimeRejected(testCase)
    folder = tempname;
    mkdir(folder);
    cleanupObj = onCleanup(@() cleanupDirectory(folder)); %#ok<NASGU>
    filePath = fullfile(folder, 'bad_time.txt');
    content = sprintf(['# Chunk-averaged data for SPID\n', ...
        '# Units microscale\n', ...
        '# Name bin1d Kind spatial\n', ...
        '# Timestep Number-of-chunks Total-count\n', ...
        '# Chunk Coord1 Ncount\n', ...
        '# Time 1oops\n', ...
        '100 1 1\n', ...
        '1 0 1\n']);
    writeTextFile(filePath, content);
    verifyError(testCase, @() read_chunk_step_fast(filePath, ...
        'SelectBy', 'Index', 'Index', 1, 'ProgressMode', 'off'), ...
        'read_chunk_step_fast:BadPhysicalTime');
end

function testChunkReadersRejectCrossLineColumnBleed(testCase)
    folder = tempname;
    mkdir(folder);
    cleanupObj = onCleanup(@() cleanupDirectory(folder)); %#ok<NASGU>
    filePath = fullfile(folder, 'bad_rows.txt');
    content = sprintf(['# header 1\n', ...
        '# header 2\n', ...
        '# Chunk Coord1 Ncount\n', ...
        '100 2 2\n', ...
        '1 0 1\n', ...
        '2 1\n', ...
        '200 1 1\n', ...
        '1 2 3\n']);
    writeTextFile(filePath, content);
    verifyError(testCase, @() read_chunk_step_fast(filePath, ...
        'SelectBy', 'Index', 'Index', 1, 'ProgressMode', 'off'), ...
        'read_chunk_step_fast:BadDataRow');
    verifyError(testCase, @() read_bin_chunk(filePath, ...
        'ProgressMode', 'off'), 'read_bin_chunk:BadDataRow');
end

function testChunkReadersAcceptBlankLinesInsideBlock(testCase)
    folder = tempname;
    mkdir(folder);
    cleanupObj = onCleanup(@() cleanupDirectory(folder)); %#ok<NASGU>
    filePath = fullfile(folder, 'blank_rows.txt');
    content = sprintf(['# header 1\n', ...
        '# header 2\n', ...
        '# Chunk Coord1 Ncount\n', ...
        '100 2 2\n', ...
        '1 0 1\n\n', ...
        '2 1 1\n']);
    writeTextFile(filePath, content);
    fast = read_chunk_step_fast(filePath, 'SelectBy', 'Index', ...
        'Index', 1, 'ProgressMode', 'off');
    full = read_bin_chunk(filePath, 'ProgressMode', 'off');
    verifyEqual(testCase, fast.data, [1 0 1; 2 1 1]);
    verifyEqual(testCase, full.steps(1).data, fast.data);
end

function testChunkReadersRejectExtraFrameSummaryColumns(testCase)
    folder = tempname;
    mkdir(folder);
    cleanupObj = onCleanup(@() cleanupDirectory(folder)); %#ok<NASGU>
    filePath = fullfile(folder, 'bad_summary.txt');
    content = sprintf(['# header 1\n', ...
        '# header 2\n', ...
        '# Chunk Coord1 Ncount\n', ...
        '100 1 1 999\n', ...
        '1 0 1\n']);
    writeTextFile(filePath, content);
    verifyError(testCase, @() read_chunk_step_fast(filePath, ...
        'SelectBy', 'Index', 'Index', 1, 'ProgressMode', 'off'), ...
        'read_chunk_step_fast:BadBlockHeader');
    verifyError(testCase, @() read_bin_chunk(filePath, ...
        'ProgressMode', 'off'), 'read_bin_chunk:BadBlockHeader');
end

function testChunk3dInputAutoResolution(testCase)
    folder = tempname;
    mkdir(folder);
    cleanupObj = onCleanup(@() cleanupDirectory(folder)); %#ok<NASGU>
    filePath = fullfile(folder, 'bin3d_case.txt');
    writeTextFile(filePath, 'placeholder');
    request = pd_create_request('chunk');
    request.baseDir = folder;
    request.chunk.dimension = '3d';
    verifyEqual(testCase, pd_resolve_input_file(request), filePath);
end

function testExplicitChunk2dRejectsOneDimensionalInput(testCase)
    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    request.chunk.dimension = '2d';
    verifyError(testCase, @() postdata_run(request), ...
        'postdata:PreflightMissingVariable');
end

function testChunkGradientAndStrainRate(testCase)
    gradient = makeRequest(testCase, 'chunk', 'bin1d_strain_rate_dx_1.txt');
    gradient.chunk.variable = 'gradient';
    gradient.chunk.gradientVariable = 'c_rho';
    gradient.chunk.coordScale = 2;
    resultGradient = postdata_run(gradient);
    verifyEqual(testCase, resultGradient.dimension, '1d');
    verifyEqual(testCase, resultGradient.x, (0:2:8).', 'AbsTol', 1e-12);
    verifyEqual(testCase, resultGradient.value, 0.5 .* ones(5, 1), ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, resultGradient.derivedSource.type, 'gradient1d');

    strain = makeRequest(testCase, 'chunk', 'bin1d_strain_rate_dx_1.txt');
    strain.chunk.variable = 'strainRate';
    strain.chunk.strainRateVelocityComponent = 'vz';
    strain.chunk.strainRateDensityVariable = 'c_rho';
    resultStrain = postdata_run(strain);
    density = (2:6).';
    velocity = (1:2:9).';
    verifyEqual(testCase, resultStrain.value, 2 + velocity ./ density, ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, resultStrain.derivedSource.type, 'strainRate1d');
end

function testClusterAnalysis(testCase)
    request = makeRequest(testCase, 'cluster', 'cluster_chunk_test.txt');
    request.analysisOptions = {'Dx', 0.5, 'Dim', 2};
    result = postdata_run(request);
    verifyEqual(testCase, result.analysisType, 'cluster');
    verifyEqual(testCase, result.selectedRows, 5);
    verifyEqual(testCase, numel(result.diameter), 5);
end

function testClusterPhysicalParticleSizeModels(testCase)
    request = makeRequest(testCase, 'cluster', 'cluster_chunk_test.txt');
    request.resultLevel = 'full';
    request.analysisOptions = {'Dim', 3, 'ParticleVolume', 0.08};
    result = postdata_run(request);
    count = result.filteredData(:, result.colIndex.Ncount);
    expected = 2 .* ((3 .* count .* 0.08) ./ (4 .* pi)).^(1 / 3);
    verifyEqual(testCase, result.diameter, expected, 'AbsTol', 1e-12);
    verifyEqual(testCase, result.sizeModel.mode, 'explicit-volume-3d');

    request.analysisOptions = {'Dim', 3, 'ParticleVolume', 0.08, ...
        'ThinDirectionThickness', 0.2};
    result = postdata_run(request);
    expected = 2 .* sqrt((count .* 0.08 ./ 0.2) ./ pi);
    verifyEqual(testCase, result.diameter, expected, 'AbsTol', 1e-12);
    verifyEqual(testCase, result.sizeModel.mode, 'quasi-2d-projection');
end

function testVelocityAnalysis(testCase)
    request = makeRequest(testCase, 'vx', 'vx_chunk_test.txt');
    result = postdata_run(request);
    verifyEqual(testCase, result.analysisType, 'vx');
    verifyEqual(testCase, numel(result.velocity), 5);
end

function testCurrentSpidFieldVelocityAndUnits(testCase)
    request = makeRequest(testCase, 'vx', 'spid_field_vx.txt');
    result = postdata_run(request);
    verifyEqual(testCase, result.velocityVar, 'Coord1');
    verifyEqual(testCase, result.velocity, [-10;0;10], 'AbsTol', 1e-12);
    totalIndex = find(strcmp(result.densityVars, ...
        'massArealDensity'), 1, 'first');
    verifyEqual(testCase, result.cumulativeDensity(:, totalIndex), ...
        [3;2.5;1.5], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.velocityFactor, 10);
    verifyEqual(testCase, result.densityFactor, 1);
    verifyTrue(testCase, result.usedFileUnitMetadata);
    verifyEqual(testCase, result.chunkKind, 'field');
end

function testCurrentSpidVelocityUsesRequestedDisplayUnit(testCase)
    request = makeRequest(testCase, 'vx', 'spid_field_vx.txt');
    request.analysisOptions = {'VelocityUnit', 'm/s'};
    result = postdata_run(request);
    verifyEqual(testCase, result.velocity, [-10000;0;10000], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.velocityFactor, 10000);
    verifyEqual(testCase, result.velocityUnit, 'm/s');
end

function testCurrentSpidClusterColumnAliases(testCase)
    request = makeRequest(testCase, 'cluster', 'spid_cluster.txt');
    request.analysisOptions = {'Dx', 0.5, 'Dim', 2, ...
        'Range_c_x', [1 2], 'XVarForMean', 'c_x'};
    result = postdata_run(request);
    verifyEqual(testCase, result.selectedRows, 2);
    verifyEqual(testCase, result.filteredData(:, result.colIndex.x), [1;2]);
    verifyFalse(testCase, isempty(result.meanByBin.centers));
    verifyEqual(testCase, result.chunkKind, 'cluster');
end

function testLegacyClusterColumnAliasesAreBidirectional(testCase)
    request = makeRequest(testCase, 'cluster', 'cluster_chunk_test.txt');
    request.analysisOptions = {'Dx', 0.5, 'Dim', 2, ...
        'XVarForMean', 'x'};
    result = postdata_run(request);
    verifyFalse(testCase, isempty(result.meanByBin.centers));
    verifyTrue(testCase, any(result.meanByBin.count > 0));
end

function testSpidUnitConversionTable(testCase)
    si = pd_spid_unit_info('si');
    cgs = pd_spid_unit_info('centimeter_gram_microsecond');
    micro = pd_spid_unit_info('microscale');
    verifyEqual(testCase, si.lengthUmPerUnit, 1e6);
    verifyEqual(testCase, si.arealDensityMgCm2PerUnit, 100);
    verifyEqual(testCase, cgs.velocityKmSPerUnit, 10);
    verifyEqual(testCase, cgs.arealDensityMgCm2PerUnit, 1000);
    verifyEqual(testCase, micro.lengthUmPerUnit, 10);
    verifyEqual(testCase, micro.arealDensityMgCm2PerUnit, 1);
end

function testMassVAliasesAndEditableCoordinate(testCase)
    request = makeRequest(testCase, 'mass-v', 'vx_chunk_test.txt');
    request.analysisOptions = {'VelocityVar', 'Chunk', 'VelocityFactor', 0.002, ...
        'VelocityLabel', 'Particle velocity', 'VelocityUnit', 'm/s', ...
        'CumulativeDirection', 'low-to-high'};
    result = postdata_run(request);
    verifyEqual(testCase, request.analysisType, 'vx');
    verifyEqual(testCase, result.velocity, (0:4).' .* 0.002, 'AbsTol', 1e-12);
    verifyEqual(testCase, result.velocityLabel, 'Particle velocity');
    verifyEqual(testCase, result.velocityUnit, 'm/s');
    verifyEqual(testCase, result.cumulativeDirection, 'low-to-high');
end

function testMassXHighToLowCumulativeDistribution(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    result = postdata_run(request);
    verifyEqual(testCase, result.analysisType, 'massx');
    verifyEqual(testCase, result.simulationType, 'SPH');
    verifyEqual(testCase, result.massModel, 'planar-2d-count');
    verifyEqual(testCase, result.coordinate, (0:5:20).', 'AbsTol', 1e-12);
    verifyEqual(testCase, result.localParticleCount, (1:5).');
    verifyEqual(testCase, result.density, (1:5).' .* 0.01, ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.cumulativeDensity, ...
        [0.15; 0.14; 0.12; 0.09; 0.05], 'AbsTol', 1e-12);
    verifyEqual(testCase, result.cumulativeDirection, 'high-to-low');
    verifyTrue(testCase, all(result.isMonotonicDecreasing));
    verifyLessThanOrEqual(testCase, ...
        max(diff(result.cumulativeDensity)), 1e-12);
end

function testCurrentSpidMassXUsesLengthMetadata(testCase)
    request = makeRequest(testCase, 'massx', 'spid_spatial1d.txt');
    result = postdata_run(request);
    verifyEqual(testCase, result.coordinate, (0:5:20).', 'AbsTol', 1e-12);
    verifyEqual(testCase, result.rawLengthUnitUm, 10);
    verifyTrue(testCase, result.usedFileUnitMetadata);
    verifyEqual(testCase, result.chunkKind, 'spatial');
    verifyTrue(testCase, all(result.isMonotonicDecreasing));
end

function testCurrentSpidMassXUsesRequestedCoordinateUnit(testCase)
    request = makeRequest(testCase, 'massx', 'spid_spatial1d.txt');
    request.analysisOptions = pd_upsert_option( ...
        request.analysisOptions, 'CoordinateUnit', 'mm');
    result = postdata_run(request);
    verifyEqual(testCase, result.coordinate, (0:0.005:0.02).', ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.coordinateFactor, 0.01, ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.coordinateUnit, 'mm');
end

function testMassXRangeAndFixedDecreasingDirection(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    request.analysisOptions = pd_upsert_option( ...
        request.analysisOptions, 'CoordinateFactor', 20);
    request.analysisOptions = pd_upsert_option( ...
        request.analysisOptions, 'CoordinateRange', [10 30]);
    result = postdata_run(request);
    verifyEqual(testCase, result.coordinate, [10;20;30], 'AbsTol', 1e-12);
    verifyEqual(testCase, result.cumulativeDensity, ...
        [0.09;0.07;0.04], 'AbsTol', 1e-12);
    verifyTrue(testCase, result.isMonotonicDecreasing);

    request.analysisOptions = {'CumulativeDirection', 'low-to-high'};
    verifyError(testCase, @() postdata_run(request), ...
        'postdata:InvalidOption');
end

function testMassXTwoDimensionalFullWidthAndSlices(testCase)
    request = makeRequest(testCase, 'massx', ...
        'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.resultLevel = 'full';
    request.analysisOptions = {'ChunkDim', '2d', 'InitialDensity', 2, ...
        'ParticleSpacing', 0.2, 'RawLengthUnitUm', 10, ...
        'CoordinateFactor', 10};
    fullResult = postdata_run(request);
    verifyEqual(testCase, fullResult.chunkDimension, '2d');
    verifyEqual(testCase, fullResult.coordinate, [0;10], 'AbsTol', 1e-12);
    verifyEqual(testCase, fullResult.localParticleCount, [2;4], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, fullResult.density, [0.08;0.16], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, fullResult.cumulativeDensity, [0.24;0.16], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, fullResult.particleMassModelValue, 0.08, ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, fullResult.particleLineMassGPerCm, 8e-8, ...
        'AbsTol', 1e-18);
    verifyEqual(testCase, fullResult.sliceDefinitions.coveredWidthRaw, 2, ...
        'AbsTol', 1e-12);

    request.analysisOptions = {'ChunkDim', '2d', 'InitialDensity', 2, ...
        'ParticleSpacing', 0.2, 'RawLengthUnitUm', 10, ...
        'CoordinateFactor', 10, 'SliceCentersY', [0 1], ...
        'SliceWidthsY', 1};
    slicedResult = postdata_run(request);
    verifyEqual(testCase, slicedResult.localParticleCount, [2 0;1 3], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, slicedResult.density, [0.16 0;0.08 0.24], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, slicedResult.cumulativeDensity, ...
        [0.24 0.24;0.08 0.24], 'AbsTol', 1e-12);
    verifyTrue(testCase, all(slicedResult.isMonotonicDecreasing));
    verifyEqual(testCase, numel(slicedResult.sliceDefinitions), 2);
end

function testMassXTwoDimensionalFractionalSliceOverlap(testCase)
    request = makeRequest(testCase, 'massx', ...
        'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.analysisOptions = {'ChunkDim', '2d', 'InitialDensity', 2, ...
        'ParticleSpacing', 0.2, 'RawLengthUnitUm', 10, ...
        'CoordinateFactor', 10, 'SliceCentersY', 0.5, ...
        'SliceWidthsY', 1};
    result = postdata_run(request);
    verifyEqual(testCase, result.localParticleCount, [1;2], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.density, [0.08;0.16], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.sliceBoundaryMethod, ...
        'fractional-bin-overlap');
end

function testMassXThreeDimensionalParticleMassAndAreaNormalization(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    request.analysisOptions = {'ChunkDim', '1d', 'SphDimension', 3, ...
        'InitialDensity', 2, 'ParticleSpacing', 0.1, ...
        'RawLengthUnitUm', 10, 'TransverseWidth', 1, ...
        'OutOfPlaneWidth', 2, 'CoordinateFactor', 10};
    result = postdata_run(request);
    verifyEqual(testCase, result.massModel, 'volumetric-3d-count');
    verifyEqual(testCase, result.particleMassModelValue, 0.002, ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.particleMassG, 2e-12, ...
        'AbsTol', 1e-22);
    verifyTrue(testCase, isnan(result.particleLineMassGPerCm));
    verifyEqual(testCase, result.density, (1:5).' .* 0.001, ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, result.cumulativeDensity, ...
        [0.015;0.014;0.012;0.009;0.005], 'AbsTol', 1e-12);
end

function testSharedMassDistributionSupportsPrefixedColumns(testCase)
    step = struct();
    step.data = [1 2 3; 4 5 9];
    step.colIndex = struct('v_mass1ArealDensity', 1, ...
        'v_mass2ArealDensity', 2, 'v_massArealDensity', 3);
    step.validVarNames = {'v_mass1ArealDensity','v_mass2ArealDensity', ...
        'v_massArealDensity'};
    distribution = pd_cumulative_areal_density(step, [0;1], {}, ...
        'high-to-low', 'test_mass');
    totalIndex = find(strcmp(distribution.densityVars, ...
        'v_massArealDensity'), 1, 'first');
    verifyEqual(testCase, distribution.density(:, totalIndex), [3;9]);
    verifyEqual(testCase, distribution.cumulativeDensity(:, totalIndex), [12;9]);
end

function testSharedMassDistributionAggregatesDuplicateCoordinates(testCase)
    step = struct();
    step.data = [1;2;3];
    step.colIndex = struct('massArealDensity', 1);
    step.validVarNames = {'massArealDensity'};
    distribution = pd_cumulative_areal_density(step, [0;0;1], {}, ...
        'high-to-low', 'test_mass');
    verifyEqual(testCase, distribution.coordinate, [0;1]);
    verifyEqual(testCase, distribution.density, [3;3]);
    verifyEqual(testCase, distribution.cumulativeDensity, [6;3]);
end

function testNegativeMassDensityCannotBreakMonotonicCumulative(testCase)
    step = struct();
    step.data = [-1;2;3];
    step.colIndex = struct('massArealDensity', 1);
    step.validVarNames = {'massArealDensity'};
    warning('off', 'test_mass:NegativeDensityClipped');
    cleanupObj = onCleanup(@() warning('on', ...
        'test_mass:NegativeDensityClipped')); %#ok<NASGU>
    distribution = pd_cumulative_areal_density(step, [0;1;2], {}, ...
        'high-to-low', 'test_mass', 'clip');
    verifyEqual(testCase, distribution.density, [0;2;3]);
    verifyEqual(testCase, distribution.cumulativeDensity, [5;5;3]);
    verifyTrue(testCase, all(distribution.isMonotonicDecreasing));
    verifyEqual(testCase, distribution.negativeValueCount, 1);
    verifyError(testCase, @() pd_cumulative_areal_density(step, [0;1;2], {}, ...
        'high-to-low', 'test_mass', 'error'), 'test_mass:NegativeDensity');
end

function testNetworkAnalysis(testCase)
    request = makeRequest(testCase, 'network2d', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original'};
    result = postdata_run(request);
    verifyEqual(testCase, result.analysisType, 'network2d');
    verifyEqual(testCase, result.global.porosity, 0.25, 'AbsTol', 1e-12);
    verifyEqual(testCase, nnz(result.validMask), 4);
end

function testNetworkRendererDisabledContract(testCase)
    request = makeRequest(testCase, 'network2d', ...
        'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.resultLevel = 'full';
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original'};
    result = postdata_run(request);
    names = fieldnames(result.plots);
    for i = 1:numel(names)
        verifyEmpty(testCase, result.plots.(names{i}));
    end
end

function testNetworkRendererCompleteSnapshot(testCase)
    previousFigures = findall(0, 'Type', 'figure');
    previousVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanupObj = onCleanup(@() cleanupNetworkFigures( ...
        previousFigures, previousVisible)); %#ok<NASGU>

    request = makeRequest(testCase, 'network2d', ...
        'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.makePlots = true;
    request.resultLevel = 'full';
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original', ...
        'PositionAxis', 'both', 'ProfileAxis', 'both'};
    result = postdata_run(request);
    expected = {'phaseFig','poreLabelFig','matrixLabelFig','connectivityFig', ...
        'poreCountFig','matrixCountFig','porePositionFigX','porePositionFigY', ...
        'matrixPositionFigX','matrixPositionFigY','profileFigX','profileFigY'};
    for i = 1:numel(expected)
        verifyTrue(testCase, ishghandle(result.plots.(expected{i})), expected{i});
    end
    verifyEmpty(testCase, result.plots.evolutionFig);
end

function testRendererUsesCallerAxes(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    result = postdata_run(request);
    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    handle = pd_render_result(ax, result);
    verifyTrue(testCase, ishghandle(handle));
    verifyEqual(testCase, ancestor(handle, 'axes'), ax);
end

function testPlotDataContractMatchesRenderedViews(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.resultLevel = 'full';
    result = postdata_run(request);

    fieldData = pd_result_plot_data(result, 'field');
    verifyEqual(testCase, fieldData.ColumnNames, {'x','y',result.variableUsed});
    verifyEqual(testCase, fieldData.Values, ...
        [result.x(:), result.y(:), result.value(:)], 'AbsTol', 1e-12);
    verifyTrue(testCase, fieldData.IsRectangularGrid);

    histogramData = pd_result_plot_data(result, 'histogram');
    verifyEqual(testCase, sum(histogramData.Values(:, 2)), ...
        nnz(isfinite(result.value)), 'AbsTol', 1e-12);

    profileData = pd_result_plot_data(result, 'profile-x');
    verifyEqual(testCase, profileData.ColumnNames{1}, 'x');
    verifyEqual(testCase, size(profileData.Values, 2), 2);

    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    imageHandle = pd_render_result(ax, result, [], 'field');
    verifyEqual(testCase, get(imageHandle, 'CData'), fieldData.GridZ, ...
        'AbsTol', 1e-12);

    massRequest = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    massResult = postdata_run(massRequest);
    massData = pd_result_plot_data(massResult, 'cumulative');
    lineHandles = pd_render_result(ax, massResult, [], 'cumulative');
    verifyEqual(testCase, get(lineHandles(1), 'XData').', ...
        massData.Values(:, 1), 'AbsTol', 1e-12);
    verifyEqual(testCase, get(lineHandles(1), 'YData').', ...
        massData.Values(:, 2), 'AbsTol', 1e-12);
end

function testPlotDataCsvAndDelimitedText(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    result = postdata_run(request);
    plotData = pd_result_plot_data(result, 'cumulative');

    text = pd_plot_data_text(plotData, char(9), [1 3]);
    lines = regexp(strtrim(text), '\n', 'split');
    verifyEqual(testCase, numel(lines), 3);
    verifyNotEmpty(testCase, strfind(lines{1}, char(9))); %#ok<STRIFCND>

    textWithoutHeader = pd_plot_data_text( ...
        plotData, char(9), [1 3], false);
    linesWithoutHeader = regexp(strtrim(textWithoutHeader), '\n', 'split');
    verifyEqual(testCase, numel(linesWithoutHeader), 2);
    firstCells = regexp(linesWithoutHeader{1}, char(9), 'split');
    verifyEqual(testCase, str2double(firstCells{1}), ...
        plotData.Values(1, 1), 'AbsTol', 1e-12);
    verifyError(testCase, @() pd_plot_data_text( ...
        plotData, char(9), [1 3], 2), ...
        'postdata:plotDataText:BadHeaderFlag');

    outputDir = tempname;
    mkdir(outputDir);
    cleanupObj = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    files = pd_export_plot_data_csv(result, fullfile(outputDir, 'case'));
    views = pd_result_plot_views(result);
    verifyEqual(testCase, numel(files), numel(views));
    verifyTrue(testCase, all(cellfun(@(p) exist(p, 'file') == 2, files)));
    content = fileread(files{1});
    verifyNotEmpty(testCase, strfind(content, plotData.ColumnNames{1})); %#ok<STRIFCND>
    verifyEqual(testCase, numel(regexp(strtrim(content), '\n', 'split')), ...
        plotData.RowCount + 1);
end

function testEveryAvailableViewProvidesPlotData(testCase)
    cases = { ...
        'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt', {}; ...
        'cluster', 'cluster_chunk_test.txt', {}; ...
        'vx', 'vx_chunk_test.txt', {}; ...
        'massx', 'bin1d_dx_0.5.txt', sampleMassXOptions(); ...
        'network2d', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt', ...
            {'ThresholdN', 1, 'ProfileAxis', 'both'}};
    for i = 1:size(cases, 1)
        request = makeRequest(testCase, cases{i, 1}, cases{i, 2});
        request.resultLevel = 'full';
        request.analysisOptions = cases{i, 3};
        result = postdata_run(request);
        views = pd_result_plot_views(result);
        for j = 1:numel(views)
            plotData = pd_result_plot_data(result, views(j).Id);
            verifyEqual(testCase, plotData.ViewId, views(j).Id);
            verifyEqual(testCase, plotData.RowCount, size(plotData.Values, 1));
            verifyEqual(testCase, plotData.ColumnCount, size(plotData.Values, 2));
            verifyEqual(testCase, numel(plotData.ColumnNames), ...
                plotData.ColumnCount);
        end
    end
end

function testPlotUpdateModeParsing(testCase)
    catalog = pd_plot_option_catalog();
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'UpdateMode', 'add');
    options = pd_plot_options_from_table(catalog, data);
    verifyEqual(testCase, options.UpdateMode, 'overlay');

    data = setTableValue(data, 'UpdateMode', 'overwrite');
    options = pd_plot_options_from_table(catalog, data);
    verifyEqual(testCase, options.UpdateMode, 'replace');

    data = setTableValue(data, 'UpdateMode', 'invalid');
    verifyError(testCase, @() pd_plot_options_from_table(catalog, data), ...
        'postdata:BadPlotUpdateMode');
end

function testOptionEditorValueFormatting(testCase)
    verifyEqual(testCase, pd_format_option_editor_value([], 'numeric'), '');
    verifyEqual(testCase, pd_format_option_editor_value('[7.0]', 'numeric'), '7');
    verifyEqual(testCase, pd_format_option_editor_value('[0.1]', 'numeric'), '0.1');
    verifyEqual(testCase, ...
        pd_format_option_editor_value('[0.6 1.2 1.8]', 'numeric'), ...
        '[0.6 1.2 1.8]');

    catalog = pd_option_catalog('massx');
    data = pd_catalog_table_data(catalog);
    densityRow = find(strcmp('InitialDensity', data(:, 1)), 1, 'first');
    spacingRow = find(strcmp('ParticleSpacing', data(:, 1)), 1, 'first');
    rangeRow = find(strcmp('CoordinateRange', data(:, 1)), 1, 'first');
    verifyEqual(testCase, data{densityRow, 2}, '');
    verifyEqual(testCase, data{spacingRow, 2}, '');
    verifyEqual(testCase, data{rangeRow, 2}, '');

    data{densityRow, 2} = '[7.0]';
    data{spacingRow, 2} = '[0.1]';
    data{rangeRow, 2} = '[0 30]';
    data = pd_normalize_catalog_table_data(catalog, data);
    verifyEqual(testCase, data{densityRow, 2}, '7');
    verifyEqual(testCase, data{spacingRow, 2}, '0.1');
    verifyEqual(testCase, data{rangeRow, 2}, '[0 30]');
end

function testPlotTypedEditorUsesPopupAndCheckbox(testCase)
    assumeTrue(testCase, exist('PostDataApp', 'class') == 8, ...
        'GUI tests are intentionally excluded from the non-GUI POST_DATA sync.');
    previousVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    app = PostDataApp();
    cleanupObj = onCleanup(@() cleanupAppFigure(app, previousVisible)); %#ok<NASGU>

    plotData = get(app.PlotTable, 'Data');
    updateRow = find(strcmp('UpdateMode', plotData(:, 1)), 1, 'first');
    app.selectPlotParameter([], struct('Indices', [updateRow 2]));
    verifyEqual(testCase, char(get(app.PlotParameterPopup, 'Visible')), 'on');
    verifyEqual(testCase, get(app.PlotParameterPopup, 'String'), ...
        {'replace'; 'overlay'});

    plotData = get(app.PlotTable, 'Data');
    gridRow = find(strcmp('ShowGrid', plotData(:, 1)), 1, 'first');
    app.selectPlotParameter([], struct('Indices', [gridRow 2]));
    verifyEqual(testCase, char(get(app.PlotParameterCheckbox, 'Visible')), 'on');

    set(app.PlotPresetPopup, 'Value', 2);
    app.applyPlotPreset([], []);
    plotData = get(app.PlotTable, 'Data');
    modeRow = find(strcmp('SeriesStyleMode', plotData(:, 1)), 1, 'first');
    paletteRow = find(strcmp('ColorPalette', plotData(:, 1)), 1, 'first');
    verifyEqual(testCase, plotData{modeRow, 2}, 'monochrome');
    verifyEqual(testCase, plotData{paletteRow, 2}, 'grayscale');
end

function testRendererReplaceAndOverlay(testCase)
    request = makeRequest(testCase, 'chunk', ...
        'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    result1 = postdata_run(request);
    result2 = result1;
    result2.timestep = result1.timestep + 1;
    result2.value = result1.value .* 2;

    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    catalog = pd_plot_option_catalog();
    options = pd_plot_options_from_table(catalog, pd_catalog_table_data(catalog));

    options.UpdateMode = 'replace';
    handle1 = pd_render_result(ax, result1, options);
    verifyTrue(testCase, ishghandle(handle1));

    options.UpdateMode = 'overlay';
    handle2 = pd_render_result(ax, result2, options);
    verifyTrue(testCase, ishghandle(handle1));
    verifyTrue(testCase, ishghandle(handle2));
    proxies = findobj(ax, 'Type', 'line');
    labels = get(proxies, 'DisplayName');
    if ischar(labels), labels = {labels}; end
    verifyEqual(testCase, numel(unique(labels)), 2);

    options.UpdateMode = 'replace';
    handle3 = pd_render_result(ax, result2, options);
    verifyFalse(testCase, ishghandle(handle1));
    verifyFalse(testCase, ishghandle(handle2));
    verifyTrue(testCase, ishghandle(handle3));
end

function testCompleteOptionCatalogs(testCase)
    types = {'chunk','cluster','vx','massx','network2d'};
    minimumCounts = [6, 19, 2, 17, 34];
    for i = 1:numel(types)
        catalog = pd_option_catalog(types{i});
        verifyGreaterThanOrEqual(testCase, numel(catalog), minimumCounts(i));
        request = pd_create_request(types{i});
        data = pd_catalog_table_data(catalog);
        if strcmp(types{i}, 'massx')
            data = setTableValue(data, 'InitialDensity', '1');
            data = setTableValue(data, 'ParticleSpacing', '0.1');
            data = setTableValue(data, 'TransverseWidth', '1');
        end
        request = pd_apply_analysis_options(request, catalog, data);
        verifyEqual(testCase, request.analysisType, types{i});
    end
end

function testCatalogMatchesEveryAnalyzerOption(testCase)
    rootDir = testCase.TestData.RootDir;
    cases = { ...
        'chunk', fullfile(rootDir, 'src', 'analysis', 'analyze_chunk_field.m'); ...
        'cluster', fullfile(rootDir, 'src', 'analysis', 'cluster_postprocess.m'); ...
        'vx', fullfile(rootDir, 'src', 'analysis', 'vx_chunk_cumulative.m'); ...
        'massx', fullfile(rootDir, 'src', 'analysis', 'mass_x_cumulative.m'); ...
        'network2d', fullfile(rootDir, 'src', 'analysis', 'network2d', 'analyze_chunk_network2d.m')};
    common = {'SelectBy','Index','TimeStep','Time','SlurmPath','SlurmModuleIndex', ...
        'ProgressMode','MakePlots','MakePlot','DoPlot','PlotOptions', ...
        'CancelCallback','ProgressCallback'};
    for i = 1:size(cases, 1)
        source = fileread(cases{i, 2});
        tokens = regexp(source, 'p\.addParameter\(''([^'']+)''', 'tokens');
        analyzerNames = cellfun(@(x) x{1}, tokens, 'UniformOutput', false);
        analyzerNames = analyzerNames(~ismember(analyzerNames, common));
        catalog = pd_option_catalog(cases{i, 1});
        catalogNames = {catalog.Name};
        verifyEqual(testCase, sort(catalogNames), sort(analyzerNames), ...
            sprintf('Catalog mismatch for %s.', cases{i, 1}));
    end
end

function testSafeOptionParser(testCase)
    verifyEqual(testCase, pd_parse_option_value('[1 2.5 -3]', 'numeric'), [1 2.5 -3]);
    verifyEqual(testCase, pd_parse_option_value('true', 'logical'), true);
    verifyEqual(testCase, pd_parse_option_value('gamma,lognormal', 'celltext'), ...
        {'gamma','lognormal'});
    verifyError(testCase, @() pd_parse_option_value('system(''bad'')', 'numeric'), ...
        'postdata:BadNumericOption');
end

function testPlotConditions(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    result = postdata_run(request);
    catalog = pd_plot_option_catalog();
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'Colormap', 'jet');
    data = setTableValue(data, 'ColorLimits', '[0 2]');
    data = setTableValue(data, 'XLim', '[0 1]');
    data = setTableValue(data, 'ShowGrid', 'false');
    options = pd_plot_options_from_table(catalog, data);
    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    pd_render_result(ax, result, options);
    verifyEqual(testCase, get(ax, 'CLim'), [0 2], 'AbsTol', 1e-12);
    verifyEqual(testCase, get(ax, 'XLim'), [0 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, char(get(ax, 'XGrid')), 'off');
end

function testPublicationStyleIsFullyEditable(testCase)
    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    result = postdata_run(request);
    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    catalog = pd_plot_option_catalog();
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'ColorPalette', 'highcontrast');
    data = setTableValue(data, 'FontName', 'Arial');
    data = setTableValue(data, 'FontSize', '13');
    data = setTableValue(data, 'LineWidth', '2.4');
    data = setTableValue(data, 'LineStyle', '--');
    data = setTableValue(data, 'MarkerSymbol', 'o');
    data = setTableValue(data, 'MarkerSize', '7');
    data = setTableValue(data, 'Title', 'Paper-ready result');
    data = setTableValue(data, 'XLabel', 'Position (cm)');
    options = pd_plot_options_from_table(catalog, data);
    lineHandle = pd_render_result(ax, result, options);
    verifyEqual(testCase, get(ax, 'FontName'), 'Arial');
    verifyEqual(testCase, get(ax, 'FontSize'), 13);
    verifyEqual(testCase, get(lineHandle, 'LineWidth'), 2.4, 'AbsTol', 1e-12);
    verifyEqual(testCase, get(lineHandle, 'LineStyle'), '--');
    verifyEqual(testCase, get(lineHandle, 'Marker'), 'o');
    verifyEqual(testCase, get(get(ax, 'Title'), 'String'), 'Paper-ready result');
    verifyEqual(testCase, get(get(ax, 'XLabel'), 'String'), 'Position (cm)');

    outputCatalog = pd_output_option_catalog(tempdir);
    outputData = pd_catalog_table_data(outputCatalog);
    outputData = setTableValue(outputData, 'FigureWidthCm', '8.6');
    outputData = setTableValue(outputData, 'FigureHeightCm', '6.2');
    output = pd_output_options_from_table(outputCatalog, outputData);
    verifyEqual(testCase, output.FigureWidthCm, 8.6, 'AbsTol', 1e-12);
    verifyEqual(testCase, output.FigureHeightCm, 6.2, 'AbsTol', 1e-12);
end

function testRectangularFieldRenderModes(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    result = postdata_run(request);
    catalog = pd_plot_option_catalog();
    options = pd_plot_options_from_table(catalog, pd_catalog_table_data(catalog));
    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);

    imageHandle = pd_render_result(ax, result, options);
    verifyEqual(testCase, get(imageHandle, 'Type'), 'image');
    verifyEmpty(testCase, findobj(ax, 'Type', 'legend'));

    options.FieldRenderMode = 'scatter';
    scatterHandle = pd_render_result(ax, result, options);
    verifyEqual(testCase, get(scatterHandle, 'Type'), 'scatter');

    options.FieldRenderMode = 'contour';
    contourHandle = pd_render_result(ax, result, options);
    verifyTrue(testCase, ishghandle(contourHandle));
    [rows, columns] = pd_subplot_grid(7);
    verifyEqual(testCase, [rows, columns], [2, 4]);
end

function testResultExport(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    result = postdata_run(request);
    fig = figure('Visible', 'off');
    cleanupFig = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    pd_render_result(ax, result);

    outputDir = tempname;
    mkdir(outputDir);
    cleanupDir = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    catalog = pd_output_option_catalog(outputDir);
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'Prefix', 'integration');
    data = setTableValue(data, 'DPI', '90');
    options = pd_output_options_from_table(catalog, data);
    files = pd_export_result(result, ax, options);
    verifyEqual(testCase, numel(files), 4);
    verifyTrue(testCase, all(cellfun(@(p) exist(p, 'file') == 2, files)));
end

function testPreflightLogAndRunManifest(testCase)
    outputDir = tempname;
    mkdir(outputDir);
    cleanupObj = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.execution.logFile = fullfile(outputDir, 'run.log');
    result = postdata_run(request);
    verifyTrue(testCase, result.preflight.passed);
    verifyEqual(testCase, result.preflight.analysisType, 'chunk');
    verifyGreaterThanOrEqual(testCase, result.run.durationSeconds, 0);
    verifyEqual(testCase, exist(result.run.logFile, 'file'), 2);
    verifyEqual(testCase, result.run.request.execution.cancelCallback, '<runtime callback>');
end

function testPreflightRejectsWrongInputType(testCase)
    request = makeRequest(testCase, 'network2d', 'vx_chunk_test.txt');
    request.analysisOptions = {'ThresholdN', 1};
    verifyError(testCase, @() postdata_run(request), 'postdata:PreflightMissingVariable');
end

function testMassXPreflightRequiresParticleCount(testCase)
    request = makeRequest(testCase, 'massx', 'vx_chunk_test.txt');
    verifyError(testCase, @() postdata_run(request), ...
        'postdata:PreflightMissingVariable');
end

function testCancellationBeforeAnalysis(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.execution.cancelCallback = @() true;
    verifyError(testCase, @() postdata_run(request), 'postdata:UserCancelled');
end

function testProgrammaticRequestUsesCatalogValidation(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    request.analysisOptions = {'CumulativeDirection', 'low-to-high'};
    verifyError(testCase, @() pd_validate_request(request), ...
        'postdata:InvalidOption');

    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    request.analysisOptions = {'UnknownOption', 1};
    verifyError(testCase, @() pd_validate_request(request), ...
        'postdata:validateRequest:UnknownAnalysisOption');

    request = makeRequest(testCase, 'cluster', 'cluster_chunk_test.txt');
    request.analysisOptions = {'Dim', 2, 'dim', 3};
    verifyError(testCase, @() pd_validate_request(request), ...
        'postdata:validateRequest:DuplicateAnalysisOption');

    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    request.selection.value = 1.5;
    verifyError(testCase, @() pd_validate_request(request), ...
        'postdata:validateRequest:BadSelectionIndex');

    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    request.plotOptions = struct('UnknownStyle', 1);
    verifyError(testCase, @() postdata_run(request), ...
        'postdata:UnknownPlotOption');
end

function testAmbiguousChunkHeadersAreRejected(testCase)
    workDir = tempname;
    mkdir(workDir);
    cleanupObj = onCleanup(@() cleanupDirectory(workDir)); %#ok<NASGU>
    filePath = fullfile(workDir, 'ambiguous.txt');
    writeTextFile(filePath, sprintf([ ...
        '# header 1\n# header 2\n# Chunk Coord1 c-a c_a\n' ...
        '0 1 1\n1 0 1 2\n']));
    request = pd_create_request('chunk');
    request.baseDir = workDir;
    request.filePath = filePath;
    request.progressMode = 'off';
    request.makePlots = false;
    verifyError(testCase, @() postdata_run(request), ...
        'postdata:preflight:AmbiguousVariableName');
end

function testCorruptIndexCacheIsRebuilt(testCase)
    workDir = tempname;
    mkdir(workDir);
    cleanupObj = onCleanup(@() cleanupDirectory(workDir)); %#ok<NASGU>
    source = fullfile(testCase.TestData.FixtureDir, 'bin1d_dx_0.5.txt');
    target = fullfile(workDir, 'bin1d_dx_0.5.txt');
    copyfile(source, target);
    cachePath = [target, '.stepidx.mat'];
    writeTextFile(cachePath, 'corrupt cache');
    warningId = 'read_chunk_step_fast:InvalidIndexCache';
    previous = warning('query', warningId);
    warning('off', warningId);
    warningCleanup = onCleanup(@() warning(previous.state, warningId)); %#ok<NASGU>
    step = read_chunk_step_fast(target, 'SelectBy', 'Index', 'Index', 1, ...
        'ProgressMode', 'off');
    verifyEqual(testCase, step.stepIndex, 1);
    verifyTrue(testCase, isfield(step, 'index'));
    saved = load(cachePath, 'idx');
    verifyTrue(testCase, isfield(saved, 'idx'));
    verifyEqual(testCase, saved.idx.formatVersion, 3);
    verifyEqual(testCase, saved.idx.timesteps, step.index.timesteps);
end

function testReaderHonorsCancellation(testCase)
    filePath = fullfile(testCase.TestData.FixtureDir, 'bin1d_dx_0.5.txt');
    verifyError(testCase, @() read_chunk_step_fast(filePath, ...
        'SelectBy', 'Index', 'Index', 1, 'ProgressMode', 'off', ...
        'CancelCallback', @() true), 'postdata:UserCancelled');
end

function testSystemDiagnosticsAndDefaultOutputIsolation(testCase)
    report = pd_system_diagnostics();
    verifyTrue(testCase, report.isCompatibleMatlab);
    verifyTrue(testCase, report.passed);
    verifyTrue(testCase, all(report.moduleAvailable));
    catalog = pd_output_option_catalog();
    directoryRow = find(strcmp('Directory', {catalog.Name}), 1, 'first');
    verifyEqual(testCase, catalog(directoryRow).Default, ...
        fullfile(testCase.TestData.RootDir, 'outputs'));
end

function testFailedExportLeavesNoPartialFiles(testCase)
    request = makeRequest(testCase, 'chunk', 'bin1d_dx_0.5.txt');
    result = postdata_run(request);
    outputDir = tempname;
    mkdir(outputDir);
    cleanupObj = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    catalog = pd_output_option_catalog(outputDir);
    options = pd_output_options_from_table(catalog, pd_catalog_table_data(catalog));
    verifyError(testCase, @() pd_export_result(result, [], options), ...
        'postdata:MissingExportAxes');
    entries = dir(outputDir);
    entries = entries(~ismember({entries.name}, {'.','..'}));
    verifyEmpty(testCase, entries);
end

function testConfigurationUpgrade(testCase)
    assumeTrue(testCase, exist('pd_upgrade_app_config', 'file') == 2, ...
        'GUI configuration tests are intentionally excluded from the non-GUI POST_DATA sync.');
    oldConfig = struct('schemaVersion', 1, 'task', 'network2d', ...
        'selectionValue', '3');
    config = pd_upgrade_app_config(oldConfig);
    verifyEqual(testCase, config.schemaVersion, 2);
    verifyEqual(testCase, config.selectionValue, '3');
    verifyEqual(testCase, size(config.computeData, 1), numel(pd_option_catalog('network2d')));
    verifyEqual(testCase, size(config.plotData, 1), numel(pd_plot_option_catalog()));
    verifyEqual(testCase, size(config.outputData, 1), numel(pd_output_option_catalog()));
    verifyTrue(testCase, config.copyPlotDataHeaders);
end

function testSemanticOptionValidation(testCase)
    catalog = pd_option_catalog('network2d');
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'PositionNumBins', '2.5');
    verifyError(testCase, @() pd_validate_catalog_values(catalog, data), ...
        'postdata:InvalidOption');
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'PlotRangeX', '[10 2]');
    verifyError(testCase, @() pd_validate_catalog_values(catalog, data), ...
        'postdata:InvalidOption');
    chunkCatalog = pd_option_catalog('chunk');
    chunkData = setTableValue(pd_catalog_table_data(chunkCatalog), 'CoordScale', '0');
    verifyError(testCase, @() pd_validate_catalog_values(chunkCatalog, chunkData), ...
        'postdata:InvalidOption');
    vxCatalog = pd_option_catalog('vx');
    vxData = setTableValue(pd_catalog_table_data(vxCatalog), 'VelocityFactor', '0');
    verifyError(testCase, @() pd_validate_catalog_values(vxCatalog, vxData), ...
        'postdata:InvalidOption');
end

function testConditionalOptionsAreNotDispatched(testCase)
    catalog = pd_option_catalog('network2d');
    data = pd_catalog_table_data(catalog);
    enabled = pd_catalog_enabled_mask(catalog, data);
    evolutionStride = find(strcmp('EvolutionStride', {catalog.Name}), 1, 'first');
    verifyFalse(testCase, enabled(evolutionStride));

    request = pd_create_request('network2d');
    request = pd_apply_analysis_options(request, catalog, data);
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'EvolutionStride', 'missing'), 'missing');

    data = setTableValue(data, 'GeometryMode', 'original');
    enabled = pd_catalog_enabled_mask(catalog, data);
    cutMethod = find(strcmp('CutCellMethod', {catalog.Name}), 1, 'first');
    verifyFalse(testCase, enabled(cutMethod));
    request = pd_create_request('network2d');
    request = pd_apply_analysis_options(request, catalog, data);
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'CutCellMethod', 'missing'), 'missing');
end

function testSummaryResultRemainsRenderable(testCase)
    request = makeRequest(testCase, 'network2d', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.resultLevel = 'summary';
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original'};
    result = postdata_run(request);
    verifyEqual(testCase, result.storage.level, 'summary');
    verifyFalse(testCase, isfield(result, 'NcountGrid'));
    verifyFalse(testCase, isfield(result, 'matrix'));
    verifyTrue(testCase, isfield(result, 'poreMask'));
    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    handle = pd_render_result(ax, result);
    verifyTrue(testCase, ishghandle(handle));
end

function testAnalysisSpecificMultiViewRendering(testCase)
    request = makeRequest(testCase, 'network2d', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.resultLevel = 'full';
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original', ...
        'ProfileAxis', 'both'};
    result = postdata_run(request);
    views = pd_result_plot_views(result);
    verifyGreaterThanOrEqual(testCase, numel(views), 5);
    verifyTrue(testCase, all(ismember({'phase','pore-label','matrix-label'}, {views.Id})));
    fig = figure('Visible', 'off');
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
    for i = 1:numel(views)
        ax = axes('Parent', fig);
        handle = pd_render_result(ax, result, [], views(i).Id);
        verifyFalse(testCase, isempty(handle), views(i).Id);
        delete(ax);
    end
end

function testAppAllViewsDashboard(testCase)
    assumeTrue(testCase, exist('PostDataApp', 'class') == 8, ...
        'GUI tests are intentionally excluded from the non-GUI POST_DATA sync.');
    previousVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    app = PostDataApp();
    cleanupObj = onCleanup(@() cleanupAppFigure(app, previousVisible)); %#ok<NASGU>
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    app.CurrentResult = postdata_run(request);
    app.updatePlotViewChoices();
    app.replotCurrentResult([], []);
    verifyEqual(testCase, numel(app.PlotAxes), 3);
    verifyTrue(testCase, all(ishghandle(app.PlotAxes)));
    verifyEqual(testCase, app.PlotViewIds, ...
        {'all','cumulative','differential','particle-count'});
    verifyEqual(testCase, numel(app.PlotFigureHandles), 3);
    verifyTrue(testCase, all(ishghandle(app.PlotFigureHandles)));
    verifyTrue(testCase, all(ishghandle(app.PlotFigureAxes)));
end

function testAppPlotDataInspectorPagination(testCase)
    assumeTrue(testCase, exist('PostDataApp', 'class') == 8, ...
        'GUI tests are intentionally excluded from the non-GUI POST_DATA sync.');
    previousVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    app = PostDataApp();
    cleanupObj = onCleanup(@() cleanupAppFigure(app, previousVisible)); %#ok<NASGU>

    count = 1201;
    app.CurrentResult = struct('analysisType', 'chunk', ...
        'x', (1:count).', 'y', [], 'value', (1001:1000 + count).', ...
        'variableUsed', 'syntheticValue', 'timestep', 25, ...
        'filePath', 'synthetic.txt');
    app.updatePlotViewChoices();

    verifyEqual(testCase, get(app.PlotDataCopyHeadersCheckbox, 'Value'), 1);
    verifyEqual(testCase, app.PlotDataViewIds, {'field','histogram'});
    verifyEqual(testCase, app.CurrentPlotData.RowCount, count);
    verifyEqual(testCase, size(get(app.PlotDataTable, 'Data'), 1), 500);
    verifyEqual(testCase, app.PlotDataPage, 1);
    app.changePlotDataPage(1);
    verifyEqual(testCase, app.PlotDataPage, 2);
    verifyEqual(testCase, size(get(app.PlotDataTable, 'Data'), 1), 500);
    app.changePlotDataPage(1);
    verifyEqual(testCase, app.PlotDataPage, 3);
    verifyEqual(testCase, size(get(app.PlotDataTable, 'Data'), 1), 201);
    app.selectPlotDataRows([], struct('Indices', [1 1; 3 2]));
    verifyEqual(testCase, app.PlotDataSelectedRows, [1001 1003]);
end

function testAppExampleDataWorkflow(testCase)
    assumeTrue(testCase, exist('PostDataApp', 'class') == 8, ...
        'GUI tests are intentionally excluded from the non-GUI POST_DATA sync.');
    previousVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    app = PostDataApp();
    cleanupObj = onCleanup(@() cleanupAppFigure(app, previousVisible)); %#ok<NASGU>

    app.loadExample('chunk1d');
    request = app.buildRequest();
    verifyEqual(testCase, request.analysisType, 'chunk');
    verifyEqual(testCase, request.chunk.dimension, '1d');
    verifyEqual(testCase, exist(request.filePath, 'file'), 2);

    app.loadExample('chunk2d');
    request = app.buildRequest();
    verifyEqual(testCase, request.chunk.dimension, '2d');

    app.loadExample('cluster');
    request = app.buildRequest();
    verifyEqual(testCase, request.analysisType, 'cluster');
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, 'Dx', NaN), ...
        0.025, 'AbsTol', 1e-12);

    app.loadExample('network2d');
    request = app.buildRequest();
    verifyEqual(testCase, request.analysisType, 'network2d');
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'ProfileAxis', ''), 'both');

    app.loadExample('massx');
    request = app.buildRequest();
    verifyEqual(testCase, request.analysisType, 'massx');
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'InitialDensity', NaN), 7.3, 'AbsTol', 1e-12);
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'CoordinateFactor', NaN), 10, 'AbsTol', 1e-12);
    computeData = get(app.ComputeTable, 'Data');
    densityRow = find(strcmp('InitialDensity', computeData(:, 1)), 1, 'first');
    spacingRow = find(strcmp('ParticleSpacing', computeData(:, 1)), 1, 'first');
    transverseRow = find(strcmp('TransverseWidth', computeData(:, 1)), 1, 'first');
    verifyEqual(testCase, computeData{densityRow, 2}, '7.3');
    verifyEqual(testCase, computeData{spacingRow, 2}, '0.005');
    verifyEqual(testCase, computeData{transverseRow, 2}, '2.46');

    app.loadExample('massx2d');
    request = app.buildRequest();
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'ChunkDim', ''), '2d');
    verifyEqual(testCase, pd_get_analysis_option(request.analysisOptions, ...
        'SliceCentersY', []), [0.6 1.2 1.8], 'AbsTol', 1e-12);

    app.loadExample('vx');
    request = app.buildRequest();
    verifyEqual(testCase, request.analysisType, 'vx');
    verifyEqual(testCase, exist(request.filePath, 'file'), 2);
    verifyGreaterThanOrEqual(testCase, ...
        numel(findall(app.Figure, 'Type', 'uimenu')), 8);
end

function testGeneratedLargeOneDimensionalTwoDimensionalAndClusterData(testCase)
    generatedDir = fullfile(testCase.TestData.RootDir, 'fixtures', 'generated');

    one = pd_create_request('chunk');
    one.filePath = fullfile(generatedDir, 'large_bin1d_dx_0.025.txt');
    one.selection.value = 2;
    one.progressMode = 'off';
    one.makePlots = false;
    result1d = postdata_run(one);
    verifyEqual(testCase, result1d.dimension, '1d');
    verifyEqual(testCase, result1d.timestep, 200);
    verifyEqual(testCase, numel(result1d.x), 401);

    two = pd_create_request('chunk');
    two.filePath = fullfile(generatedDir, ...
        'large_bin2d_dx_0.05_dy_0.06_Lz_1.txt');
    two.selection.value = 3;
    two.progressMode = 'off';
    two.makePlots = false;
    result2d = postdata_run(two);
    verifyEqual(testCase, result2d.dimension, '2d');
    verifyEqual(testCase, result2d.timestep, 300);
    verifyEqual(testCase, numel(result2d.x), 2501);

    cluster = pd_create_request('cluster');
    cluster.filePath = fullfile(generatedDir, 'large_cluster_chunk.txt');
    cluster.selection.value = 2;
    cluster.progressMode = 'off';
    cluster.makePlots = false;
    cluster.analysisOptions = {'Dim', 2, 'Dx', 0.025};
    clusterResult = postdata_run(cluster);
    verifyEqual(testCase, clusterResult.timestep, 200);
    verifyEqual(testCase, clusterResult.selectedRows, 1200);
    verifyGreaterThan(testCase, max(clusterResult.diameter), ...
        min(clusterResult.diameter));

    massX = pd_create_request('mass-x');
    massX.filePath = fullfile(generatedDir, 'large_bin1d_dx_0.025.txt');
    massX.selection.value = 2;
    massX.progressMode = 'off';
    massX.makePlots = false;
    massX.analysisOptions = {'ChunkDim', '1d', 'InitialDensity', 7.3, ...
        'ParticleSpacing', 0.005, 'RawLengthUnitUm', 10, ...
        'TransverseWidth', 2.46, 'CoordinateFactor', 10};
    massXResult = postdata_run(massX);
    verifyEqual(testCase, numel(massXResult.coordinate), 401);
    verifyGreaterThan(testCase, massXResult.cumulativeDensity(1), 0);
    verifyTrue(testCase, massXResult.isMonotonicDecreasing);

    massX2d = pd_create_request('mass-x');
    massX2d.filePath = fullfile(generatedDir, ...
        'large_bin2d_dx_0.05_dy_0.06_Lz_1.txt');
    massX2d.selection.value = 2;
    massX2d.progressMode = 'off';
    massX2d.makePlots = false;
    massX2d.analysisOptions = {'ChunkDim', '2d', 'InitialDensity', 7.3, ...
        'ParticleSpacing', 0.01, 'RawLengthUnitUm', 10, ...
        'CoordinateFactor', 10, 'SliceCentersY', [0.6 1.2 1.8], ...
        'SliceWidthsY', 0.6};
    massX2dResult = postdata_run(massX2d);
    verifyEqual(testCase, size(massX2dResult.cumulativeDensity), [61 3]);
    verifyTrue(testCase, all(massX2dResult.isMonotonicDecreasing));

    massV = pd_create_request('mass-v');
    massV.filePath = fullfile(generatedDir, 'large_mass_v.txt');
    massV.selection.value = 3;
    massV.progressMode = 'off';
    massV.makePlots = false;
    massVResult = postdata_run(massV);
    verifyEqual(testCase, massVResult.timestep, 300);
    verifyEqual(testCase, numel(massVResult.velocity), 501);
    verifyGreaterThan(testCase, massVResult.cumulativeDensity(1, end), 0);
end

function testDetailCsvDispatch(testCase)
    outputDir = tempname;
    mkdir(outputDir);
    cleanupObj = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>

    chunk = struct('analysisType', 'chunk', 'x', [1;2], 'y', [], 'value', [3;4]);
    cluster = struct('analysisType', 'cluster', 'diameter', [1;2]);
    vx = struct('analysisType', 'vx', 'velocity', [1;2], 'density', [3;4], ...
        'cumulativeDensity', [7;4], 'densityVars', {{'massArealDensity'}});
    massx = struct('analysisType', 'massx', 'coordinate', [1;2], ...
        'density', [3;4], 'cumulativeDensity', [7;4], ...
        'densityVars', {{'massArealDensity'}});
    network = struct('analysisType', 'network2d', 'xCenters', [0 1], ...
        'yCenters', [0 1], 'validMask', true(2), 'poreMask', logical([1 0;0 1]));
    inputs = {chunk, cluster, vx, massx, network};
    for i = 1:numel(inputs)
        files = pd_export_detail_csv(inputs{i}, fullfile(outputDir, ['case', num2str(i)]));
        verifyGreaterThanOrEqual(testCase, numel(files), 1);
        verifyTrue(testCase, all(cellfun(@(p) exist(p, 'file') == 2, files)));
    end
end

function testChunk3dDetailCsvIncludesZ(testCase)
    request = makeRequest(testCase, 'chunk', 'spid_spatial3d.txt');
    request.chunk.variable = 'rho';
    result = postdata_run(request);
    outputDir = tempname;
    mkdir(outputDir);
    cleanupObj = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    files = pd_export_detail_csv(result, fullfile(outputDir, 'field3d'));
    verifyEqual(testCase, numel(files), 1);
    [header, values] = readNumericCsv(files{1}, 4);
    verifyEqual(testCase, header, 'x,y,z,rho');
    verifyEqual(testCase, numel(values), 4);
    verifyEqual(testCase, values{3}, result.coordZ);
end

function testResultExportIncludesAllPlotDataCsvFiles(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    result = postdata_run(request);
    outputDir = tempname;
    mkdir(outputDir);
    cleanupObj = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    catalog = pd_output_option_catalog(outputDir);
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'SaveMAT', 'false');
    data = setTableValue(data, 'SaveCSV', 'false');
    data = setTableValue(data, 'SaveDetailCSV', 'false');
    data = setTableValue(data, 'SavePlotDataCSV', 'true');
    data = setTableValue(data, 'SavePNG', 'false');
    data = setTableValue(data, 'SaveManifest', 'false');
    options = pd_output_options_from_table(catalog, data);
    files = pd_export_result(result, [], options);
    verifyEqual(testCase, numel(files), numel(pd_result_plot_views(result)));
    verifyTrue(testCase, all(cellfun(@(p) exist(p, 'file') == 2, files)));
end

function testFigPdfAndManifestExport(testCase)
    request = makeRequest(testCase, 'chunk', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    result = postdata_run(request);
    fig = figure('Visible', 'off');
    cleanupFig = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes('Parent', fig);
    pd_render_result(ax, result);
    outputDir = tempname;
    mkdir(outputDir);
    cleanupDir = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    catalog = pd_output_option_catalog(outputDir);
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'SaveMAT', 'false');
    data = setTableValue(data, 'SaveCSV', 'false');
    data = setTableValue(data, 'SavePNG', 'true');
    data = setTableValue(data, 'SaveFIG', 'true');
    data = setTableValue(data, 'SavePDF', 'true');
    data = setTableValue(data, 'SaveManifest', 'true');
    options = pd_output_options_from_table(catalog, data);
    files = pd_export_result(result, ax, options);
    verifyEqual(testCase, numel(files), 4);
    extensions = cellfun(@(p) lower(fileExtension(p)), files, 'UniformOutput', false);
    verifyTrue(testCase, all(ismember({'.png','.fig','.pdf','.txt'}, extensions)));
    manifest = files{find(strcmp('.txt', extensions), 1, 'first')};
    content = fileread(manifest);
    verifyNotEmpty(testCase, strfind(content, 'analysisType=chunk')); %#ok<STRIFCND>
    verifyNotEmpty(testCase, strfind(content, '[Request]')); %#ok<STRIFCND>
end

function testDashboardAxesExport(testCase)
    request = makeRequest(testCase, 'massx', 'bin1d_dx_0.5.txt');
    result = postdata_run(request);
    views = pd_result_plot_views(result);
    fig = figure('Visible', 'off');
    cleanupFig = onCleanup(@() close(fig)); %#ok<NASGU>
    axesList = gobjects(1, numel(views));
    for i = 1:numel(views)
        axesList(i) = axes('Parent', fig);
        pd_render_result(axesList(i), result, [], views(i).Id);
    end
    outputDir = tempname;
    mkdir(outputDir);
    cleanupDir = onCleanup(@() cleanupDirectory(outputDir)); %#ok<NASGU>
    catalog = pd_output_option_catalog(outputDir);
    data = pd_catalog_table_data(catalog);
    data = setTableValue(data, 'SaveMAT', 'false');
    data = setTableValue(data, 'SaveCSV', 'false');
    data = setTableValue(data, 'SaveDetailCSV', 'false');
    data = setTableValue(data, 'SavePNG', 'true');
    data = setTableValue(data, 'SaveManifest', 'false');
    options = pd_output_options_from_table(catalog, data);
    files = pd_export_result(result, axesList, options);
    verifyEqual(testCase, numel(files), 1);
    verifyEqual(testCase, exist(files{1}, 'file'), 2);
end

function testNetworkSpacingModule(testCase)
    options = struct('Dx', [], 'Dy', [], 'CoordScale', 1);
    filePath = fullfile(testCase.TestData.FixtureDir, 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    [dx, dy, source] = pd_network_resolve_spacing(options, filePath, [0 0.5], [0 0.5]);
    verifyEqual(testCase, [dx dy], [0.5 0.5], 'AbsTol', 1e-12);
    verifyEqual(testCase, source.dx, 'filename');
    options.Dx = 2;
    options.Dy = 3;
    options.CoordScale = 10;
    [dx, dy, source] = pd_network_resolve_spacing(options, 'plain.txt', [0 1], [0 1]);
    verifyEqual(testCase, [dx dy], [20 30]);
    verifyEqual(testCase, source.dx, 'option');
end

function testSharedStatisticsModules(testCase)
    values = [1 2 3 4];
    stats = pd_stats_summary(values, 1, 0);
    verifyEqual(testCase, stats.mean, 2.5, 'AbsTol', 1e-12);
    verifyEqual(testCase, stats.quantileValue, [1.15 1.75 2.5 3.25 3.85], 'AbsTol', 1e-12);
    verifyEqual(testCase, pd_stats_moment_ratio(values, 2, 1), 3, 'AbsTol', 1e-12);
    emptyStats = pd_stats_summary([]);
    verifyTrue(testCase, isnan(emptyStats.mean));
end

function testEvolutionModuleIntegration(testCase)
    request = makeRequest(testCase, 'network2d', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original', ...
        'EnableEvolution', true, 'EvolutionSelectBy', 'Index', ...
        'EvolutionRange', [1 2], 'EvolutionStride', 1};
    result = postdata_run(request);
    evolution = result.stats.evolution;
    verifyEqual(testCase, evolution.axis, 'index');
    verifyEqual(testCase, evolution.stepIndex, [1 2]);
    verifyEqual(testCase, numel(evolution.geometry.phi), 2);
    verifyTrue(testCase, isnan(evolution.transition.deltaPhi(1)));
end

function testDigitalTopologyOpenCases(testCase)
    fullMask = true(3);
    [labels, ~, ~] = pd_network_label_components(fullMask, true(3), 'open');
    topology = pd_network_compute_topology(fullMask, 'open', max(labels(:)), false);
    verifyEqual(testCase, [topology.beta0 topology.beta1 topology.chi], [1 0 1]);

    ring = true(3);
    ring(2, 2) = false;
    [labels, ~, ~] = pd_network_label_components(ring, true(3), 'open');
    topology = pd_network_compute_topology(ring, 'open', max(labels(:)), false);
    verifyEqual(testCase, [topology.beta0 topology.beta1 topology.chi], [1 1 0]);
end

function testDigitalTopologyPeriodicCases(testCase)
    cylinder = true(1, 4);
    [labels, wrapsX, wrapsY] = pd_network_label_components(cylinder, true(size(cylinder)), 'periodic-x');
    topology = pd_network_compute_topology(cylinder, 'periodic-x', max(labels(:)), false);
    verifyEqual(testCase, topology.beta1, 1);
    verifyTrue(testCase, wrapsX(1));
    verifyFalse(testCase, wrapsY(1));

    torus = true(2, 2);
    [labels, wrapsX, wrapsY] = pd_network_label_components(torus, true(size(torus)), 'periodic-xy');
    topology = pd_network_compute_topology(torus, 'periodic-xy', max(labels(:)), false);
    verifyEqual(testCase, [topology.beta0 topology.beta1 topology.beta2 topology.chi], [1 2 1 0]);
    verifyTrue(testCase, wrapsX(1));
    verifyTrue(testCase, wrapsY(1));
end

function testPeriodicSeamConnectivityWithoutWinding(testCase)
    mask = logical([1 0 1]);
    [openLabels, ~, ~] = pd_network_label_components(mask, true(size(mask)), 'open');
    verifyEqual(testCase, max(openLabels(:)), 2);
    [periodicLabels, wrapsX, wrapsY] = ...
        pd_network_label_components(mask, true(size(mask)), 'periodic-x');
    verifyEqual(testCase, max(periodicLabels(:)), 1);
    verifyFalse(testCase, wrapsX(1));
    verifyFalse(testCase, wrapsY(1));
end

function testMorphologyDisabledAndUnavailable(testCase)
    unavailable = struct('bwskel', false, 'bwmorph', false, ...
        'distanceAvailable', false, 'skeletonAvailable', false, ...
        'imageProcessingAvailable', false);
    phaseMask = logical([1 1; 0 0]);
    disabled = pd_network_build_morphology_stats( ...
        phaseMask, ~phaseMask, 2, 2, 1, 1, 'open', false, unavailable);
    enabled = pd_network_build_morphology_stats( ...
        phaseMask, ~phaseMask, 2, 2, 1, 1, 'open', true, unavailable);

    verifyFalse(testCase, disabled.network.pore.enabled);
    verifyEqual(testCase, disabled.network.pore.note, ...
        'Skeleton graph is disabled (EnableSkeletonGraph=false).');
    verifyTrue(testCase, isnan(disabled.network.pore.skeletonLength));
    verifyTrue(testCase, enabled.network.pore.enabled);
    verifyFalse(testCase, enabled.network.pore.available);
    verifyEqual(testCase, enabled.network.pore.note, ...
        'Skeleton analysis requires Image Processing Toolbox (bwskel or bwmorph).');
    verifyEqual(testCase, enabled.thickness.matrix.note, ...
        'Thickness analysis requires Image Processing Toolbox (bwdist plus bwskel/bwmorph).');
end

function testMorphologyEmptyPhaseContract(testCase)
    available = struct('bwskel', true, 'bwmorph', true, ...
        'distanceAvailable', true, 'skeletonAvailable', true, ...
        'imageProcessingAvailable', true);
    morphology = pd_network_build_morphology_stats( ...
        false(2, 3), false(2, 3), 0, 0, 2, 3, ...
        'periodic-xy', true, available);

    verifyEqual(testCase, morphology.network.pore.skeletonLength, 0);
    verifyEqual(testCase, morphology.network.pore.branchPoints, 0);
    verifyEqual(testCase, morphology.network.pore.endPoints, 0);
    verifyEqual(testCase, morphology.network.pore.branchDensity, 0);
    verifyEqual(testCase, morphology.network.pore.note, 'Phase is empty.');
    verifyEqual(testCase, morphology.thickness.matrix.sampleCount, 0);
    verifyEqual(testCase, morphology.thickness.matrix.note, 'Matrix phase is empty.');
end

function testMorphologyAnisotropicSkeletonLength(testCase)
    toolbox = imageToolboxAvailability();
    assumeTrue(testCase, toolbox.skeletonAvailable);
    horizontal = false(3, 4);
    horizontal(2, :) = true;
    vertical = false(4, 3);
    vertical(:, 2) = true;

    hStats = pd_network_build_morphology_stats( ...
        horizontal, false(size(horizontal)), 24, 0, 2, 3, ...
        'open', true, toolbox);
    vStats = pd_network_build_morphology_stats( ...
        vertical, false(size(vertical)), 24, 0, 2, 3, ...
        'open', true, toolbox);

    verifyEqual(testCase, hStats.network.pore.skeletonLength, 6, 'AbsTol', 1e-12);
    verifyEqual(testCase, vStats.network.pore.skeletonLength, 9, 'AbsTol', 1e-12);
    verifyEqual(testCase, hStats.network.pore.endPoints, 2);
    verifyEqual(testCase, vStats.network.pore.endPoints, 2);
end

function testComponentSizeFilterContract(testCase)
    [diameter, area, mask] = pd_network_filter_component_sizes( ...
        [NaN 1 2 3 0], [9 10 20 30 40], [1.5 3]);
    verifyEqual(testCase, mask, [false; false; true; true; false]);
    verifyEqual(testCase, diameter, [2; 3]);
    verifyEqual(testCase, area, [20; 30]);
end

function testSharedDistributionEmptyBinModes(testCase)
    values = [1 1 3 3];
    zeroMode = pd_distribution_build(values, [1 4], 1, {}, 'zero');
    nanMode = pd_distribution_build(values, [1 4], 1, {}, 'nan');
    removeMode = pd_distribution_build(values, [1 4], 1, {}, 'remove');

    verifyEqual(testCase, zeroMode.edges, [1 2 3 4], 'AbsTol', 1e-12);
    verifyEqual(testCase, zeroMode.count, [2 0 2]);
    verifyEqual(testCase, zeroMode.probability, [0.5 0 0.5], 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(zeroMode.rawPdf .* diff(zeroMode.edges)), ...
        1, 'AbsTol', 1e-12);
    verifyEqual(testCase, zeroMode.cdfX, [1; 1; 3; 3]);
    verifyEqual(testCase, zeroMode.cdfProbability, [0.25; 0.5; 0.75; 1], ...
        'AbsTol', 1e-12);
    verifyTrue(testCase, isnan(nanMode.count(2)));
    verifyTrue(testCase, isnan(nanMode.pdf(2)));
    verifyEqual(testCase, removeMode.centers, [1.5 3.5], 'AbsTol', 1e-12);
    verifyEqual(testCase, removeMode.count, [2 2]);
    verifyEqual(testCase, removeMode.keptBins, [true false true]);
end

function testSharedDistributionFits(testCase)
    data = pd_distribution_build([1 2 3 4 5], [], 1, ...
        {'powerlaw', 'gamma', 'lognormal'}, 'zero');

    verifyTrue(testCase, data.fit.powerlaw.enabled);
    verifyTrue(testCase, data.fit.gamma.enabled);
    verifyTrue(testCase, data.fit.lognormal.enabled);
    verifyGreaterThan(testCase, data.fit.powerlaw.params.alpha, 1);
    verifyGreaterThan(testCase, data.fit.gamma.params.shape, 0);
    verifyGreaterThan(testCase, data.fit.gamma.params.scale, 0);
    verifyGreaterThan(testCase, data.fit.lognormal.params.sigma, 0);
    verifyEqual(testCase, data.fit.gamma.count, ...
        data.fit.gamma.pdf .* data.totalCount .* data.binSize, 'AbsTol', 1e-12);
    verifyEqual(testCase, data.fit.lognormal.probability, ...
        data.fit.lognormal.pdf .* data.binSize, 'AbsTol', 1e-12);
end

function testSharedDistributionOptionNormalization(testCase)
    verifyEqual(testCase, pd_distribution_normalize_empty_bin_mode('NA'), 'nan');
    verifyEqual(testCase, pd_distribution_normalize_empty_bin_mode('drop'), 'remove');
    verifyEqual(testCase, pd_distribution_normalize_fit_types( ...
        'power-law; gam; ln; gamma', 'FitTypes'), ...
        {'powerlaw', 'gamma', 'lognormal'});
    verifyEmpty(testCase, pd_distribution_normalize_fit_types('off', 'FitTypes'));
    verifyError(testCase, @() pd_distribution_normalize_fit_types( ...
        'unknown', 'FitTypes', 'test:BadFitTypes'), 'test:BadFitTypes');
end

function testDirectionalProfilesBinaryBins(testCase)
    poreMask = logical([1 0; 1 0]);
    validMask = true(2, 2);
    [labelGrid, wrapX, wrapY] = ...
        pd_network_label_components(poreMask, validMask, 'open');
    [comp, ~, ~] = pd_network_compute_component_stats( ...
        labelGrid, poreMask, validMask, [0.5 1.5], [0.5 1.5], ...
        1, 1, 'open', wrapX, wrapY, []);
    phase = struct();
    phase.components = comp;
    phase.connectivity = pd_network_build_directional_connectivity(comp, 'open', 2);
    opt = profileOptions('both', 2);

    profiles = pd_network_build_directional_profiles( ...
        poreMask, validMask, [0.5 1.5], [0.5 1.5], 1, 1, ...
        opt, phase, []);

    verifyEqual(testCase, profiles.x.validArea, [2 2], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.x.poreArea, [2 0], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.x.matrixArea, [0 2], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.x.porosity, [1 0], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.x.interfaceLength, [0 2], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.x.connectivityFlag, [true false]);
    verifyEqual(testCase, profiles.y.poreArea, [1 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.y.matrixArea, [1 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.y.porosity, [0.5 0.5], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.y.interfaceLength, [1 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, profiles.y.connectivityFlag, [false false]);
end

function testDirectionalProfilesCutCellConservation(testCase)
    ncount = [0 0.8; 0.2 2];
    validMask = true(2, 2);
    poreMask = ncount < 1;
    geomOpt = geometryOptions('cutcell');
    cutCell = pd_network_build_geometry(ncount, validMask, poreMask, ...
        ~poreMask, [0.5 1.5], [0.5 1.5], 1, 1, geomOpt);
    [labelGrid, wrapX, wrapY] = ...
        pd_network_label_components(poreMask, validMask, 'open');
    [comp, ~, ~] = pd_network_compute_component_stats( ...
        labelGrid, poreMask, validMask, [0.5 1.5], [0.5 1.5], ...
        1, 1, 'open', wrapX, wrapY, cutCell.pore);
    phase = struct();
    phase.components = comp;
    phase.connectivity = pd_network_build_directional_connectivity( ...
        comp, 'open', sum(cutCell.pore.areaGrid(:)));
    profiles = pd_network_build_directional_profiles( ...
        poreMask, validMask, [0.5 1.5], [0.5 1.5], 1, 1, ...
        profileOptions('both', 2), phase, cutCell);

    poreArea = sum(cutCell.pore.areaGrid(:));
    matrixArea = sum(cutCell.matrix.areaGrid(:));
    verifyEqual(testCase, sum(profiles.x.poreArea), poreArea, 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(profiles.y.poreArea), poreArea, 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(profiles.x.matrixArea), matrixArea, 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(profiles.y.matrixArea), matrixArea, 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(profiles.x.interfaceLength), ...
        cutCell.interfaceLength, 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(profiles.y.interfaceLength), ...
        cutCell.interfaceLength, 'AbsTol', 1e-12);
end

function testDirectionalConnectivityPeriodicCriterion(testCase)
    phaseMask = logical([1 1 1; 0 0 0]);
    validMask = true(2, 3);
    [labelGrid, wrapX, wrapY] = ...
        pd_network_label_components(phaseMask, validMask, 'periodic-x');
    [comp, ~, ~] = pd_network_compute_component_stats( ...
        labelGrid, phaseMask, validMask, 1:3, 1:2, 1, 1, ...
        'periodic-x', wrapX, wrapY, []);
    connectivity = pd_network_build_directional_connectivity( ...
        comp, 'periodic-x', 3);

    verifyEqual(testCase, connectivity.x.criterion, 'wraps-periodic-boundary');
    verifyTrue(testCase, connectivity.x.isConnected);
    verifyEqual(testCase, connectivity.x.componentIds, 1);
    verifyEqual(testCase, connectivity.x.totalConnectedArea, 3, 'AbsTol', 1e-12);
    verifyEqual(testCase, connectivity.x.connectedAreaFraction, 1, 'AbsTol', 1e-12);
    verifyEqual(testCase, connectivity.y.criterion, 'touches-both-open-boundaries');
    verifyFalse(testCase, connectivity.y.isConnected);
end

function testComponentStatsOpenBoundary(testCase)
    phaseMask = logical([1 1 0; 0 0 1]);
    validMask = true(2, 3);
    [labelGrid, wrapX, wrapY] = ...
        pd_network_label_components(phaseMask, validMask, 'open');
    [comp, rank, status] = pd_network_compute_component_stats( ...
        labelGrid, phaseMask, validMask, [1 3 5], [1.5 4.5], 2, 3, ...
        'open', wrapX, wrapY, []);

    verifyEqual(testCase, comp.cellCount, [2; 1]);
    verifyEqual(testCase, comp.area, [12; 6], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.centroidX, [2; 5], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.centroidY, [1.5; 4.5], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.bboxWidth, [4; 2], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.bboxHeight, [3; 3], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.perimeterOpen, [16; 10], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.interfacePerimeter, [8; 5], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.shapeFactor, ...
        [4 * pi * 12 / 16^2; 4 * pi * 6 / 10^2], 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.touchesLeft, [true; false]);
    verifyEqual(testCase, comp.touchesRight, [false; true]);
    verifyEqual(testCase, rank.label, [1; 2]);
    verifyEqual(testCase, rank.rankByLabel, [1; 2]);
    verifyFalse(testCase, status.percolatesX);
    verifyFalse(testCase, status.percolatesY);
    verifyFalse(testCase, status.wrapsX);
    verifyFalse(testCase, status.wrapsY);
end

function testComponentStatsPeriodicStatus(testCase)
    phaseMask = logical([1 1 1; 0 0 0]);
    validMask = true(2, 3);
    [labelGrid, wrapX, wrapY] = ...
        pd_network_label_components(phaseMask, validMask, 'periodic-x');
    [comp, ~, status] = pd_network_compute_component_stats( ...
        labelGrid, phaseMask, validMask, 1:3, 1:2, 1, 1, ...
        'periodic-x', wrapX, wrapY, []);

    verifyTrue(testCase, isnan(status.percolatesX));
    verifyTrue(testCase, status.wrapsX);
    verifyFalse(testCase, status.percolatesY);
    verifyFalse(testCase, status.wrapsY);
    verifyEqual(testCase, comp.perimeterOpen, 6, 'AbsTol', 1e-12);
    verifyEqual(testCase, comp.interfacePerimeter, 3, 'AbsTol', 1e-12);
end

function testEmptyComponentStatsContract(testCase)
    phaseMask = false(2, 2);
    [comp, rank, status] = pd_network_compute_component_stats( ...
        zeros(2, 2), phaseMask, true(2, 2), 1:2, 1:2, 1, 1, ...
        'open', false(0, 1), false(0, 1), []);

    verifyEmpty(testCase, comp.label);
    verifyEqual(testCase, comp.ownerGrid, zeros(2, 2));
    verifyEmpty(testCase, rank.label);
    verifyFalse(testCase, status.percolatesX);
    verifyFalse(testCase, status.percolatesY);
end

function testBinaryGeometryModule(testCase)
    ncount = [0 2; 2 2];
    valid = true(2);
    pore = ncount < 1;
    options = geometryOptions('original');
    geometry = pd_network_build_geometry(ncount, valid, pore, ~pore, ...
        [0.5 1.5], [0.5 1.5], 1, 1, options);
    verifyEqual(testCase, geometry.pore.fraction, double(pore));
    verifyEqual(testCase, geometry.matrix.fraction, double(~pore));
    verifyEqual(testCase, geometry.interfaceLength, 2, 'AbsTol', 1e-12);
end

function testCutCellGeometryConservation(testCase)
    ncount = [0 0.8; 0.2 2];
    valid = true(2);
    pore = ncount < 1;
    options = geometryOptions('cutcell');
    geometry = pd_network_build_geometry(ncount, valid, pore, ~pore, ...
        [0.5 1.5], [0.5 1.5], 1, 1, options);
    total = geometry.pore.fraction(valid) + geometry.matrix.fraction(valid);
    verifyEqual(testCase, total, ones(size(total)), 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, geometry.interfaceLength, 0);
    partial = geometry.pore.fraction(valid);
    verifyTrue(testCase, any(partial > 0 & partial < 1));
    verifyEqual(testCase, geometry.fallbackCellCount, 0);
end

function testCutCellFallbackModes(testCase)
    ncount = NaN;
    valid = true;
    pore = false;
    options = geometryOptions('cutcell');
    options.CutCellFallback = 'binary';
    geometry = pd_network_build_geometry(ncount, valid, pore, ~pore, 0.5, 0.5, 1, 1, options);
    verifyEqual(testCase, geometry.fallbackCellCount, 1);
    verifyEqual(testCase, geometry.matrix.fraction, 1);
    options.CutCellFallback = 'error';
    verifyError(testCase, @() pd_network_build_geometry( ...
        ncount, valid, pore, ~pore, 0.5, 0.5, 1, 1, options), ...
        'analyze_chunk_network2d:CutCellFallbackRequired');
end

function testCutCellFixtureReference(testCase)
    request = makeRequest(testCase, 'network2d', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
    request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'cutcell'};
    result = postdata_run(request);
    verifyEqual(testCase, result.global.interfaceLength, 0.87610275633, 'AbsTol', 1e-10);
    verifyEqual(testCase, result.cutCell.fallbackCellCount, 0);
    verifyEqual(testCase, result.global.poreArea + result.global.matrixArea, ...
        result.global.validArea, 'AbsTol', 1e-12);
end

function request = makeRequest(testCase, analysisType, fileName)
    request = pd_create_request(analysisType);
    request.baseDir = testCase.TestData.FixtureDir;
    request.filePath = fileName;
    request.progressMode = 'off';
    request.makePlots = false;
    if strcmp(pd_normalize_analysis_type(analysisType), 'massx')
        request.analysisOptions = sampleMassXOptions();
    end
end

function options = sampleMassXOptions()
    options = {'ChunkDim', '1d', 'InitialDensity', 1, ...
        'ParticleSpacing', 0.1, 'RawLengthUnitUm', 10, ...
        'TransverseWidth', 1, 'CoordinateFactor', 10};
end

function data = setTableValue(data, name, value)
    row = find(strcmp(name, data(:, 1)), 1, 'first');
    if isempty(row)
        error('test:MissingOption', 'Missing option %s.', name);
    end
    data{row, 2} = value;
end

function cleanupDirectory(pathValue)
    if exist(pathValue, 'dir')
        rmdir(pathValue, 's');
    end
end

function writeTextFile(filePath, content)
    fid = fopen(filePath, 'w');
    if fid < 0
        error('test:FileOpenFailed', 'Cannot create test file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', content);
end

function [header, values] = readNumericCsv(filePath, columnCount)
    fid = fopen(filePath, 'r');
    if fid < 0
        error('test:FileOpenFailed', 'Cannot read test file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    header = fgetl(fid);
    format = repmat('%f', 1, columnCount);
    values = textscan(fid, format, 'Delimiter', ',');
end

function cleanupNetworkFigures(previousFigures, previousVisible)
    set(0, 'DefaultFigureVisible', previousVisible);
    currentFigures = findall(0, 'Type', 'figure');
    for i = 1:numel(currentFigures)
        keep = false;
        for j = 1:numel(previousFigures)
            if isequal(currentFigures(i), previousFigures(j))
                keep = true;
                break;
            end
        end
        if ~keep && ishghandle(currentFigures(i))
            close(currentFigures(i));
        end
    end
end

function cleanupAppFigure(app, previousVisible)
    set(0, 'DefaultFigureVisible', previousVisible);
    if ~isempty(app) && ishghandle(app.Figure)
        close(app.Figure);
    end
end

function extension = fileExtension(pathValue)
    [~, ~, extension] = fileparts(pathValue);
end

function options = geometryOptions(mode)
    options = struct('GeometryMode', mode, 'CutCellMethod', 'plic', ...
        'CutCellFallback', 'binary', 'CutCellPlotRefinement', 1, ...
        'ThresholdN', 1, 'Boundary', 'open');
end

function options = profileOptions(axisName, nBins)
    options = struct('ProfileAxis', axisName, 'ProfileNumBins', nBins, ...
        'ProfileRangeX', [], 'ProfileRangeY', []);
end

function toolbox = imageToolboxAvailability()
    toolbox = struct();
    toolbox.bwskel = (exist('bwskel', 'file') == 2);
    toolbox.bwmorph = (exist('bwmorph', 'file') == 2);
    toolbox.distanceAvailable = (exist('bwdist', 'file') == 2);
    toolbox.skeletonAvailable = toolbox.bwskel || toolbox.bwmorph;
    toolbox.imageProcessingAvailable = ...
        toolbox.distanceAvailable && toolbox.skeletonAvailable;
end
