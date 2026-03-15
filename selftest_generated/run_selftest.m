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

set(0,'DefaultFigureVisible','off');
addpath(fileparts(baseDir));

fprintf('BASE|%s\n', baseDir); fprintf(logFid, 'BASE|%s\n', baseDir);

slurmPath = get_slurm_txt_fullpath(baseDir);
fprintf('SLURM|%s\n', slurmPath); fprintf(logFid, 'SLURM|%s\n', slurmPath);

chunkPath = fullfile(baseDir, 'bin2d_dx_0.5_dy_0.5_Lz_1.txt');
step = read_chunk_step_fast(chunkPath, 'SelectBy', 'Time', 'Time', 10, 'SlurmPath', slurmPath);
idxRho = step.colIndex.c_rho;
idxMass = step.colIndex.v_mass;
rhoSum = sum(step.data(isfinite(step.data(:,idxRho)), idxRho));
fprintf('CHUNK_STEP|timestep=%g|stepIndex=%d|rows=%d\n', step.timestep, step.stepIndex, size(step.data,1));
fprintf(logFid, 'CHUNK_STEP|timestep=%g|stepIndex=%d|rows=%d\n', step.timestep, step.stepIndex, size(step.data,1));
fprintf('CHUNK_RHO_SUM|%.6f\n', rhoSum); fprintf(logFid, 'CHUNK_RHO_SUM|%.6f\n', rhoSum);
fprintf('CHUNK_ROW3_MASS_ISNAN|%d\n', isnan(step.data(3, idxMass))); fprintf(logFid, 'CHUNK_ROW3_MASS_ISNAN|%d\n', isnan(step.data(3, idxMass)));

outChunkRaw = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'Variable', 'c_rho', 'SelectBy', 'Time', 'Time', 20, ...
    'SlurmPath', slurmPath, 'DoPlot', false);
fprintf('RUN_CHUNK_RAW|timestep=%g|n=%d\n', outChunkRaw.timestep, numel(outChunkRaw.value));
fprintf(logFid, 'RUN_CHUNK_RAW|timestep=%g|n=%d\n', outChunkRaw.timestep, numel(outChunkRaw.value));

dV = 0.5 * 0.5 * 1.0;
outChunkVM = run_analysis('chunk', 'BaseDir', baseDir, 'ChunkDim', '2d', ...
    'Variable', 'vonMisesS', 'SelectBy', 'TimeStep', 'TimeStep', 200, ...
    'dV', dV, 'DoPlot', false);
vm = outChunkVM.value;
vmFin = vm(isfinite(vm));
fprintf('RUN_CHUNK_VM|timestep=%g|min=%.6f|max=%.6f\n', outChunkVM.timestep, min(vmFin), max(vmFin));
fprintf(logFid, 'RUN_CHUNK_VM|timestep=%g|min=%.6f|max=%.6f\n', outChunkVM.timestep, min(vmFin), max(vmFin));

outCluster = run_analysis('cluster', 'BaseDir', baseDir, 'ClusterFile', 'cluster_chunk_test.txt', ...
    'SelectBy', 'Time', 'Time', 10, 'SlurmPath', slurmPath, ...
    'ClusterOptions', {'MakePlots', false, 'Dx', 0.5, 'Dim', 2});
fprintf('RUN_CLUSTER|timestep=%g|selectedRows=%d|diameterMean=%.6f\n', ...
    outCluster.timestep, outCluster.selectedRows, outCluster.stats.mean);
fprintf(logFid, 'RUN_CLUSTER|timestep=%g|selectedRows=%d|diameterMean=%.6f\n', ...
    outCluster.timestep, outCluster.selectedRows, outCluster.stats.mean);

outVx = run_analysis('vx', 'BaseDir', baseDir, 'VxFile', 'vx_chunk_test.txt', ...
    'SelectBy', 'Time', 'Time', 20, 'SlurmPath', slurmPath, ...
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
