% run_selftest.m
baseDir = pwd;
outputLog = fullfile(baseDir, 'selftest_output.txt');
if exist(outputLog, 'file')
    delete(outputLog);
end
logFid = fopen(outputLog, 'w');
if logFid < 0
    error('selftest:LogOpenFail', 'Cannot open selftest_output.txt');
end

oldFigureVisible = get(0, 'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');
figureVisibleCleanup = onCleanup(@() restoreFigureVisible(oldFigureVisible)); %#ok<NASGU>
addpath(fileparts(baseDir));

fprintf('BASE|%s\n', baseDir); fprintf(logFid, 'BASE|%s\n', baseDir);

slurmPath = get_slurm_txt_fullpath(baseDir);
fprintf('SLURM|%s\n', slurmPath); fprintf(logFid, 'SLURM|%s\n', slurmPath);

expectError(@() run_analysis('cluster', 'BaseDir', baseDir, 'ClusterFile', 'cluster_chunk_test', ...
    'ClusterOptions', {'MakePlots', false, 'Dx', 0.5, 'Dim', 2}), ...
    'run_analysis:MissingTxtSuffix');
expectError(@() run_analysis('vx', 'BaseDir', baseDir, 'VxFile', 'vx_chunk_test', ...
    'VxOptions', {'MakePlot', false}), ...
    'run_analysis:MissingTxtSuffix');
expectError(@() run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'ChunkFile', 'bin2d_dx_0.5_dy_0.5_Lz_1', 'Variable', 'c_rho', 'DoPlot', false), ...
    'run_analysis:MissingTxtSuffix');
fprintf('EXPECT_BAD_SUFFIX|cluster=1|vx=1|chunk=1\n');
fprintf(logFid, 'EXPECT_BAD_SUFFIX|cluster=1|vx=1|chunk=1\n');

chunkPath = fullfile(baseDir, 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
progressTempDir = tempname;
mkdir(progressTempDir);
progressCleanup = onCleanup(@() cleanupTempDir(progressTempDir)); %#ok<NASGU>
progressChunkPath = fullfile(progressTempDir, 'progress_chunk_test.txt');
copyOk = copyfile(chunkPath, progressChunkPath);
assertTrue(copyOk == 1, 'Failed to prepare temp chunk file for progress tests.');

progressText1 = evalc('read_chunk_step_fast(progressChunkPath, ''Index'', 1, ''ProgressMode'', ''console'');');
assertContainsText(progressText1, 'Building index', ...
    'First indexed read should print building-index progress.');
assertContainsText(progressText1, '100%', ...
    'First indexed read should report 100%% completion.');
progressText2 = evalc('read_chunk_step_fast(progressChunkPath, ''Index'', 1, ''ProgressMode'', ''console'');');
assertTrue(isempty(strtrim(progressText2)), ...
    'Cached read_chunk_step_fast call should not print progress.');
progressText3 = evalc('read_bin_chunk(progressChunkPath, ''ProgressMode'', ''console'');');
assertContainsText(progressText3, 'Reading file', ...
    'read_bin_chunk should print reading progress.');
assertContainsText(progressText3, '100%', ...
    'read_bin_chunk should report 100%% completion.');
progressText4 = evalc('read_bin_chunk(progressChunkPath, ''ProgressMode'', ''console'');');
assertContainsText(progressText4, 'Reading file', ...
    'Repeated read_bin_chunk call should still print progress.');
progressTextAuto = evalc('read_bin_chunk(progressChunkPath, ''ProgressMode'', ''auto'');');
if usejava('desktop') && usejava('awt')
    assertTrue(isempty(strtrim(progressTextAuto)) || ~isempty(strfind(progressTextAuto, 'Reading file')), ...
        'ProgressMode=auto should use waitbar or console output on desktop MATLAB.');
else
    assertContainsText(progressTextAuto, 'Reading file', ...
        'ProgressMode=auto should fall back to console output in batch mode.');
end
progressTextOff = evalc('read_bin_chunk(progressChunkPath, ''ProgressMode'', ''off'');');
assertTrue(isempty(strtrim(progressTextOff)), ...
    'read_bin_chunk with ProgressMode=off should not print progress.');
runOffText = evalc(['run_analysis(''chunk'', ''BaseDir'', baseDir, ''ChunkDim'', ''2d'', ' ...
    '''ChunkFile'', ''bin2d_dx_0.5_dy_0.5_Lz_1.txt'', ''Variable'', ''c_rho'', ' ...
    '''Index'', 1, ''ProgressMode'', ''off'', ''DoPlot'', false);']);
assertTrue(isempty(strtrim(runOffText)), ...
    'run_analysis with ProgressMode=off should not print progress.');
fprintf('PROGRESS_TESTS_OK\n');
fprintf(logFid, 'PROGRESS_TESTS_OK\n');

step = read_chunk_step_fast(chunkPath, 'SelectBy', 'Time', 'Time', 10, ...
    'SlurmPath', slurmPath, 'ProgressMode', 'off');
idxRho = step.colIndex.c_rho;
idxMass = step.colIndex.v_mass;
rhoSum = sum(step.data(isfinite(step.data(:,idxRho)), idxRho));
fprintf('CHUNK_STEP|timestep=%g|stepIndex=%d|rows=%d\n', step.timestep, step.stepIndex, size(step.data,1));
fprintf(logFid, 'CHUNK_STEP|timestep=%g|stepIndex=%d|rows=%d\n', step.timestep, step.stepIndex, size(step.data,1));
fprintf('CHUNK_RHO_SUM|%.6f\n', rhoSum); fprintf(logFid, 'CHUNK_RHO_SUM|%.6f\n', rhoSum);
fprintf('CHUNK_ROW3_MASS_ISNAN|%d\n', isnan(step.data(3, idxMass))); fprintf(logFid, 'CHUNK_ROW3_MASS_ISNAN|%d\n', isnan(step.data(3, idxMass)));

outChunkRaw = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'Variable', 'c_rho', 'SelectBy', 'Time', 'Time', 20, ...
    'SlurmPath', slurmPath, 'ProgressMode', 'off', 'DoPlot', false);
fprintf('RUN_CHUNK_RAW|timestep=%g|n=%d\n', outChunkRaw.timestep, numel(outChunkRaw.value));
fprintf(logFid, 'RUN_CHUNK_RAW|timestep=%g|n=%d\n', outChunkRaw.timestep, numel(outChunkRaw.value));

dV = 0.5 * 0.5 * 1.0;
outChunkVM = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'Variable', 'vonMisesS', 'SelectBy', 'TimeStep', 'TimeStep', 200, ...
    'dV', dV, 'ProgressMode', 'off', 'DoPlot', false);
vm = outChunkVM.value;
vmFin = vm(isfinite(vm));
fprintf('RUN_CHUNK_VM|timestep=%g|min=%.6f|max=%.6f\n', outChunkVM.timestep, min(vmFin), max(vmFin));
fprintf(logFid, 'RUN_CHUNK_VM|timestep=%g|min=%.6f|max=%.6f\n', outChunkVM.timestep, min(vmFin), max(vmFin));

refChunk = compute_temp_stress_chunk(chunkPath, 'SelectBy', 'Time', 'Time', 10, ...
    'SlurmPath', slurmPath, 'ProgressMode', 'off', 'dV', dV);
expectedDensity = step.data(:, idxMass) ./ dV .* 1.66053906660;
outChunkVelocity = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'ChunkFile', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt', ...
    'Variable', 'velocity', 'SelectBy', 'Time', 'Time', 10, ...
    'SlurmPath', slurmPath, 'ProgressMode', 'off', 'DoPlot', false);
outChunkPressure = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'ChunkFile', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt', ...
    'Variable', 'pressure', 'SelectBy', 'Time', 'Time', 10, ...
    'SlurmPath', slurmPath, 'ProgressMode', 'off', 'dV', dV, 'DoPlot', false);
outChunkDensity = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'ChunkFile', 'bin2d_dx_0.5_dy_0.5_Lz_1.txt', ...
    'Variable', 'density', 'SelectBy', 'Time', 'Time', 10, ...
    'SlurmPath', slurmPath, 'ProgressMode', 'off', 'dV', dV, 'DoPlot', false);
assertArrayApprox(outChunkVelocity.value, refChunk.speed, 1e-12, ...
    'Derived velocity should match compute_temp_stress_chunk speed.');
assertArrayApprox(outChunkPressure.value, refChunk.pressure, 1e-12, ...
    'Derived pressure should match compute_temp_stress_chunk pressure.');
assertArrayApprox(refChunk.density, expectedDensity, 1e-12, ...
    'compute_temp_stress_chunk density should divide every bin by the scalar dV.');
assertArrayApprox(outChunkDensity.value, expectedDensity, 1e-12, ...
    'Derived density should match the manual mass/volume conversion.');
assertTrue(sum(isfinite(outChunkDensity.value)) == sum(isfinite(expectedDensity)), ...
    'Density should stay finite for every non-empty bin when dV is scalar.');
fprintf('RUN_CHUNK_DERIVED|velocityMax=%.6f|pressureMean=%.6f|densityMean=%.6f\n', ...
    max(outChunkVelocity.value(isfinite(outChunkVelocity.value))), ...
    mean(outChunkPressure.value(isfinite(outChunkPressure.value))), ...
    mean(outChunkDensity.value(isfinite(outChunkDensity.value))));
fprintf(logFid, 'RUN_CHUNK_DERIVED|velocityMax=%.6f|pressureMean=%.6f|densityMean=%.6f\n', ...
    max(outChunkVelocity.value(isfinite(outChunkVelocity.value))), ...
    mean(outChunkPressure.value(isfinite(outChunkPressure.value))), ...
    mean(outChunkDensity.value(isfinite(outChunkDensity.value))));

outCluster = run_analysis('cluster', 'BaseDir', baseDir, 'ClusterFile', 'cluster_chunk_test.txt', ...
    'SelectBy', 'Time', 'Time', 10, 'SlurmPath', slurmPath, ...
    'ProgressMode', 'off', ...
    'ClusterOptions', {'MakePlots', false, 'Dx', 0.5, 'Dim', 2});
assertApprox(outCluster.stats.mean, mean(outCluster.diameter), 1e-12, ...
    'Default cluster mean should remain arithmetic mean.');
fprintf('RUN_CLUSTER|timestep=%g|selectedRows=%d|diameterMean=%.6f\n', ...
    outCluster.timestep, outCluster.selectedRows, outCluster.stats.mean);
fprintf(logFid, 'RUN_CLUSTER|timestep=%g|selectedRows=%d|diameterMean=%.6f\n', ...
    outCluster.timestep, outCluster.selectedRows, outCluster.stats.mean);

rangeDia = [0.95 1.35];
outClusterFiltered = run_analysis('cluster', 'BaseDir', baseDir, ...
    'ClusterFile', 'cluster_chunk_test.txt', ...
    'SelectBy', 'Time', 'Time', 10, 'SlurmPath', slurmPath, ...
    'ProgressMode', 'off', ...
    'ClusterOptions', {'MakePlots', true, 'Dx', 0.5, 'Dim', 2, ...
                       'Range_diameter', rangeDia, ...
                       'MeanPowerM', 2, 'MeanPowerN', 1, ...
                       'MeanNumBins', 2, 'HistNumBins', 4});
expectedSelected = sum(outCluster.diameter >= rangeDia(1) & outCluster.diameter <= rangeDia(2));
assertTrue(outClusterFiltered.selectedRows == expectedSelected, ...
    'Range_diameter should reduce selected rows consistently.');
expectedCustomMean = sum(outClusterFiltered.diameter .^ 2) / sum(outClusterFiltered.diameter);
assertApprox(outClusterFiltered.stats.mean, expectedCustomMean, 1e-12, ...
    'Custom cluster mean should use sum(d^m)/sum(d^n).');
assertTrue(outClusterFiltered.meanDefinition.m == 2 && outClusterFiltered.meanDefinition.n == 1, ...
    'meanDefinition should record the configured powers.');

xVals = outClusterFiltered.filteredData(:, outClusterFiltered.colIndex.c_x);
expectedBinMeans = computeBinnedMomentMeans(xVals, outClusterFiltered.diameter, ...
    outClusterFiltered.meanByBin.edges, 2, 1);
validBins = ~isnan(expectedBinMeans);
assertApprox(max(abs(outClusterFiltered.meanByBin.meanDiameter(validBins) - expectedBinMeans(validBins))), ...
    0, 1e-12, 'Binned diameter mean should use the same power-ratio definition.');

assertTrue(ishandle(outClusterFiltered.plots.countFig), 'countFig should be returned when MakePlots=true.');
assertTrue(ishandle(outClusterFiltered.plots.probFig), 'probFig should be returned when MakePlots=true.');
assertTrue(ishandle(outClusterFiltered.plots.cdfFig), 'cdfFig should be returned when MakePlots=true.');
assertTrue(ishandle(outClusterFiltered.plots.meanFig), 'meanFig should be returned when MakePlots=true.');
assertTrue(outClusterFiltered.plots.histFig == outClusterFiltered.plots.countFig, ...
    'histFig should remain as a compatibility alias of countFig.');
fprintf('RUN_CLUSTER_FILTERED|timestep=%g|selectedRows=%d|diameterMean=%.6f\n', ...
    outClusterFiltered.timestep, outClusterFiltered.selectedRows, outClusterFiltered.stats.mean);
fprintf(logFid, 'RUN_CLUSTER_FILTERED|timestep=%g|selectedRows=%d|diameterMean=%.6f\n', ...
    outClusterFiltered.timestep, outClusterFiltered.selectedRows, outClusterFiltered.stats.mean);
fprintf('RUN_CLUSTER_PLOTS|countFig=%d|probFig=%d|histAlias=%d\n', ...
    ishandle(outClusterFiltered.plots.countFig), ...
    ishandle(outClusterFiltered.plots.probFig), ...
    outClusterFiltered.plots.histFig == outClusterFiltered.plots.countFig);
fprintf(logFid, 'RUN_CLUSTER_PLOTS|countFig=%d|probFig=%d|histAlias=%d\n', ...
    ishandle(outClusterFiltered.plots.countFig), ...
    ishandle(outClusterFiltered.plots.probFig), ...
    outClusterFiltered.plots.histFig == outClusterFiltered.plots.countFig);
close([outClusterFiltered.plots.countFig, outClusterFiltered.plots.probFig, ...
       outClusterFiltered.plots.cdfFig, outClusterFiltered.plots.meanFig]);

outVx = run_analysis('vx', 'BaseDir', baseDir, 'VxFile', 'vx_chunk_test.txt', ...
    'SelectBy', 'Time', 'Time', 20, 'SlurmPath', slurmPath, ...
    'ProgressMode', 'off', ...
    'VxOptions', {'MakePlot', false, 'VelocityFactor', 0.001});
fprintf('RUN_VX|timestep=%g|nVelocity=%d|nPlotted=%d\n', ...
    outVx.timestep, numel(outVx.velocity), numel(outVx.plottedVars));
fprintf(logFid, 'RUN_VX|timestep=%g|nVelocity=%d|nPlotted=%d\n', ...
    outVx.timestep, numel(outVx.velocity), numel(outVx.plottedVars));
if ~isempty(outVx.plottedVars)
    fprintf('RUN_VX_PLOTTED|%s\n', strjoin(outVx.plottedVars, ','));
    fprintf(logFid, 'RUN_VX_PLOTTED|%s\n', strjoin(outVx.plottedVars, ','));
else
    fprintf('RUN_VX_PLOTTED|<none>\n');
    fprintf(logFid, 'RUN_VX_PLOTTED|<none>\n');
end

fprintf('SELFTEST_OK\n');
fprintf(logFid, 'SELFTEST_OK\n');
fclose(logFid);

function assertTrue(cond, msg)
if ~cond
    error('selftest:AssertionFailed', '%s', msg);
end
end

function assertApprox(actual, expected, tol, msg)
if nargin < 4
    msg = 'Values are not within tolerance.';
end
if ~(isfinite(actual) && isfinite(expected) && abs(actual - expected) <= tol)
    error('selftest:AssertionFailed', '%s (actual=%g, expected=%g, tol=%g)', ...
        msg, actual, expected, tol);
end
end

function assertArrayApprox(actual, expected, tol, msg)
if nargin < 4
    msg = 'Arrays are not within tolerance.';
end
if ~isequal(size(actual), size(expected))
    error('selftest:AssertionFailed', '%s (size mismatch).', msg);
end
sameNaN = isequal(isnan(actual), isnan(expected));
if ~sameNaN
    error('selftest:AssertionFailed', '%s (NaN mask mismatch).', msg);
end
mask = isfinite(actual) & isfinite(expected);
if any(mask(:))
    maxErr = max(abs(actual(mask) - expected(mask)));
else
    maxErr = 0;
end
if maxErr > tol
    error('selftest:AssertionFailed', '%s (maxErr=%g, tol=%g)', msg, maxErr, tol);
end
end

function assertContainsText(txt, pattern, msg)
if isempty(strfind(txt, pattern))
    error('selftest:AssertionFailed', '%s (missing pattern: %s)', msg, pattern);
end
end

function expectError(fun, expectedId)
try
    fun();
catch ME
    if strcmp(ME.identifier, expectedId)
        return;
    end
    error('selftest:UnexpectedError', ...
        'Expected error %s but got %s (%s).', expectedId, ME.identifier, ME.message);
end
error('selftest:MissingExpectedError', 'Expected error %s was not raised.', expectedId);
end

function meanVals = computeBinnedMomentMeans(x, d, edges, meanPowerM, meanPowerN)
binId = discretize(x, edges);
meanVals = nan(1, numel(edges)-1);
for i = 1:numel(meanVals)
    idx = (binId == i);
    if any(idx)
        meanVals(i) = sum(d(idx) .^ meanPowerM) / sum(d(idx) .^ meanPowerN);
    end
end
end

function cleanupTempDir(dirPath)
if exist(dirPath, 'dir')
    try
        rmdir(dirPath, 's');
    catch
        % ignore cleanup failures
    end
end
end

function restoreFigureVisible(oldValue)
try
    set(0, 'DefaultFigureVisible', oldValue);
catch
    % ignore restore failures
end
end
