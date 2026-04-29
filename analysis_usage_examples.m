% analysis_usage_examples.m
% =============================================================
% Unified usage guide for run_analysis(taskType, ...)
% =============================================================
%
% taskType = 'chunk' | 'cluster' | 'vx' | 'network2d'
%
% -------------------------------------------------------------
% Common parameters
% -------------------------------------------------------------
% 'BaseDir'          : case directory (default: pwd)
% 'SelectBy'         : 'Index' | 'TimeStep' | 'Time'
% 'Index'            : block index when SelectBy='Index'
% 'TimeStep'         : timestep when SelectBy='TimeStep'
% 'Time'             : physical time when SelectBy='Time'
% 'SlurmPath'        : slurm file path (optional; auto from BaseDir for Time mode)
% 'SlurmModuleIndex' : module index in slurm Step/CPU blocks (default: 1)
% 'ProgressMode'     : 'auto' | 'off' | 'waitbar' | 'console'
%
% -------------------------------------------------------------
% Chunk-specific additions
% -------------------------------------------------------------
% 'ChunkDim'         : '1d' or '2d'
% 'ChunkFile'        : optional explicit .txt file path/name
% 'Variable'         : raw column or derived variable
% 'dV'               : required for pressure/density/stress-related derived vars
% 'CoordScale'       : multiply chunk coordinates by this factor (default 1)
% 'CoordRangeX/Y'    : optional scaled coordinate ranges for chunk outputs
% 'DoPlot'           : true/false
% 'PlotOptions'      : cell, forwarded to plot_line1d/plot_cloud2d
%
% -------------------------------------------------------------
% Cluster-specific additions
% -------------------------------------------------------------
% 'ClusterFile'      : optional explicit .txt file path/name
% 'ClusterOptions'   : cell, forwarded to cluster_postprocess
% Note: cluster continues to use Range_c_x / Range_c_y / Range_c_z and Dx.
%
% -------------------------------------------------------------
% Vx-specific additions
% -------------------------------------------------------------
% 'VxFile'           : optional explicit .txt file path/name
% 'VxOptions'        : cell, forwarded to vx_chunk_cumulative
%
% -------------------------------------------------------------
% Network2d-specific additions
% -------------------------------------------------------------
% 'ChunkFile'        : optional explicit 2D chunk file
% 'CoordScale'       : multiply x/y coordinates and inferred Dx/Dy by this factor
% 'NetworkOptions'   : cell, forwarded to analyze_chunk_network2d
%   Required:
%   - 'ThresholdN'   : pore criterion, pore if Ncount < ThresholdN
%   Useful optional:
%   - 'Boundary'     : 'open' | 'periodic-x' | 'periodic-y' | 'periodic-xy'
%   - 'MeanPowerM/N' : generalized mean definition sum(d^m)/sum(d^n)
%   - 'PositionAxis' : axis for mean-size-vs-position distributions
%   - 'PositionRangeX/Y' : scaled centroid-coordinate ranges for position distributions
%   - 'PlotRangeX/Y' : scaled coordinate ranges used only to crop 2D network figures
%   - 'ProfileAxis'  : axis for porosity/interface/connectivity profiles
%   - 'ProfileRangeX/Y' : scaled coordinate ranges for directional profiles
%   - 'ProfileNumBins'  : number of bins for the directional profiles
%   - 'EnableSkeletonGraph' : compute skeleton metrics when toolbox is available
%   - 'EnableEvolution' : additionally scan a timestep/index/time range
%   - 'EvolutionSelectBy' : 'Index' | 'TimeStep' | 'Time'
%   - 'EvolutionRange' : [min max] range on the selected evolution axis
%   - 'EvolutionStride' : keep every Nth step in the evolution scan
%   - 'MakePlots'    : true/false
%
% Note:
%   Under periodic boundaries, "merged across the seam" and "truly connected
%   along that direction" are not identical. connectivity.x / connectivity.y
%   use the wrapped-domain winding test, while topology.beta1 is computed on
%   the wrapped 4-neighbor digital cell complex and therefore includes both
%   ordinary enclosed holes and periodic wrapping loops. Profile connectivity
%   is assigned by the centroid of each globally connected pore component.
%
% =============================================================
% Quick setup
% =============================================================
baseDir = 'C:\Users\Administrator\Desktop\WORKSPACE\codex\POST_DATA\selftest_generated';
chunkBaseDir = baseDir;
clusterBaseDir = baseDir;
vxBaseDir = baseDir;
networkBaseDir = baseDir;

chunkDim = '2d';          % '1d' or '2d'
chunkFile = '';           % leave empty to auto-match bin1d*.txt / bin2d*.txt
clusterFile = '';
vxFile = '';
networkChunkFile = '';

selectBy = 'Index';       % 'Index' | 'TimeStep' | 'Time'
index = 1;
timeStep = 200;
time = 20;
slurmPath = '';           % leave empty to auto-resolve when selectBy='Time'
slurmModuleIndex = 1;
progressMode = 'auto';
selectorArgs = {'SelectBy', selectBy, ...
                'Index', index, ...
                'TimeStep', timeStep, ...
                'Time', time, ...
                'SlurmPath', slurmPath, ...
                'SlurmModuleIndex', slurmModuleIndex, ...
                'ProgressMode', progressMode};

coordScale = 1;
networkPlotRangeX = [];
networkPlotRangeY = [];
dVMode = 'auto';          % 'auto' | 'manual'
manualDV = 4.05 * 16.1443894417461 * 4.05;

runChunkRaw = false;
runChunkDerived = false;
runCluster = false;
runVx = false;
runNetwork2d = true;

% =============================================================
% 1) Chunk analysis examples
% =============================================================
if runChunkRaw
    switch lower(strtrim(chunkDim))
        case '1d'
            out_chunk_raw = run_analysis('chunk', ...
                'BaseDir', chunkBaseDir, ...
                'ChunkDim', '1d', ...
                'ChunkFile', chunkFile, ...
                'Variable', 'Ncount', ...
                selectorArgs{:}, ...
                'CoordScale', coordScale, ...
                'CoordRangeX', [], ...
                'DoPlot', true, ...
                'PlotOptions', {'LineWidth', 1.5});

        case '2d'
            out_chunk_raw = run_analysis('chunk', ...
                'BaseDir', chunkBaseDir, ...
                'ChunkDim', '2d', ...
                'ChunkFile', chunkFile, ...
                'Variable', 'c_rho', ...
                selectorArgs{:}, ...
                'CoordScale', coordScale, ...
                'CoordRangeX', [], ...
                'CoordRangeY', [], ...
                'DoPlot', true, ...
                'PlotOptions', {'NaNMode', 'transparent', 'SmoothLevel', 0});

        otherwise
            error('analysis_usage_examples:BadChunkDim', ...
                'chunkDim must be ''1d'' or ''2d''.');
    end
end

if runChunkDerived
    if strcmpi(dVMode, 'auto')
        dV = infer_dV_from_bin_filename(chunkBaseDir, chunkDim);
    else
        dV = manualDV;
    end

    switch lower(strtrim(chunkDim))
        case '1d'
            out_chunk_derived = run_analysis('chunk', ...
                'BaseDir', chunkBaseDir, ...
                'ChunkDim', '1d', ...
                'ChunkFile', chunkFile, ...
                'Variable', 'velocity', ...
                selectorArgs{:}, ...
                'CoordScale', coordScale, ...
                'CoordRangeX', [], ...
                'dV', dV, ...
                'DoPlot', true);

        case '2d'
            out_chunk_derived = run_analysis('chunk', ...
                'BaseDir', chunkBaseDir, ...
                'ChunkDim', '2d', ...
                'ChunkFile', chunkFile, ...
                'Variable', 'vonMisesS', ...
                selectorArgs{:}, ...
                'CoordScale', coordScale, ...
                'CoordRangeX', [], ...
                'CoordRangeY', [], ...
                'dV', dV, ...
                'DoPlot', true, ...
                'PlotOptions', {'NaNMode', 'white', 'SmoothLevel', 1.0});

        otherwise
            error('analysis_usage_examples:BadChunkDim', ...
                'chunkDim must be ''1d'' or ''2d''.');
    end
end

% =============================================================
% 2) Cluster analysis example
% =============================================================
if runCluster
    out_cluster = run_analysis('cluster', ...
        'BaseDir', clusterBaseDir, ...
        'ClusterFile', clusterFile, ...
        selectorArgs{:}, ...
        'ClusterOptions', {'Dim', 2, 'Dx', 0.5, ...
                           'Range_c_x', [], ...
                           'Range_diameter', [0.1 10.0], ...
                           'MeanPowerM', 1, 'MeanPowerN', 0, ...
                           'MakePlots', true});
end

% =============================================================
% 3) Cumulative areal density (vx_chunk) example
% =============================================================
if runVx
    out_vx = run_analysis('vx', ...
        'BaseDir', vxBaseDir, ...
        'VxFile', vxFile, ...
        selectorArgs{:}, ...
        'VxOptions', {'VelocityFactor', 0.001, 'MakePlot', true});
end

% =============================================================
% 4) 2D pore/matrix network example
% =============================================================
if runNetwork2d
    out_network2d = run_analysis('network2d', ...
        'BaseDir', networkBaseDir, ...
        'ChunkDim', '2d', ...
        'ChunkFile', networkChunkFile, ...
        selectorArgs{:}, ...
        'CoordScale', coordScale, ...
        'NetworkOptions', {'ThresholdN', 1, ...
                           'Boundary', 'open', ...
                           'MeanPowerM', 1, ...
                           'MeanPowerN', 0, ...
                           'PositionAxis', 'both', ...
                           'PositionRangeX', [], ...
                           'PositionRangeY', [], ...
                           'PlotRangeX', networkPlotRangeX, ...
                           'PlotRangeY', networkPlotRangeY, ...
                           'ProfileAxis', 'x', ...
                           'ProfileRangeX', [], ...
                           'ProfileRangeY', [], ...
                           'ProfileNumBins', 24, ...
                           'EnableSkeletonGraph', true, ...
                           'EnableEvolution', false, ...
                           'EvolutionSelectBy', selectBy, ...
                           'EvolutionRange', [], ...
                           'EvolutionStride', 1, ...
                           'MakePlots', true});

    % Paper-ready snapshot metrics:
    % out_network2d.stats.topology.pore.beta0 / beta1 / chi
    % out_network2d.stats.connectivity.pore.largestFraction / percolatesX / percolatesY
    % out_network2d.stats.geometry.phi / interfaceLength / specificInterface
    % out_network2d.stats.size.pore.area / radius / diameter / hist
    % out_network2d.stats.network.pore.skeletonLength / branchPoints / endPoints
    % out_network2d.stats.thickness.matrix.min / mean / p1 / p5
    % out_network2d.stats.fragmentation.matrix.count / largestFraction
    %
    % Optional evolution scan:
    % out_network2d.stats.evolution.geometry.phi
    % out_network2d.stats.evolution.topology.pore.beta0
    % out_network2d.stats.evolution.topology.pore.beta1
    % out_network2d.stats.evolution.connectivity.pore.largestFraction
end

% =============================================================
% Helper: auto dV from bin filename metadata
% =============================================================
function dV = infer_dV_from_bin_filename(baseDir, chunkDim)
% Auto infer dV from bin filename metadata:
% use any 3 recognized keys among dx/dy/dz/Lx/Ly/Lz.

    dimTag = lower(strtrim(chunkDim));
    if strcmp(dimTag, '1d')
        pattern = 'bin1d*.txt';
    elseif strcmp(dimTag, '2d')
        pattern = 'bin2d*.txt';
    else
        error('analysis_usage_examples:BadChunkDim', 'chunkDim must be ''1d'' or ''2d''.');
    end

    files = dir(fullfile(baseDir, pattern));
    files = files(~[files.isdir]);
    if isempty(files)
        error('analysis_usage_examples:NoBinFile', ...
            'No file matches %s in %s', pattern, baseDir);
    end
    if numel(files) > 1
        names = {files.name};
        error('analysis_usage_examples:MultipleBinFiles', ...
            'Multiple files match %s in %s: %s', pattern, baseDir, strjoin(names, ', '));
    end

    fname = files(1).name;
    base = regexprep(fname, '\.txt$', '');
    parts = strsplit(base, '_');

    if numel(parts) < 3
        error('analysis_usage_examples:BadBinName', ...
            'Filename has no key-value metadata: %s', fname);
    end
    kv = parts(2:end);
    if mod(numel(kv), 2) ~= 0
        error('analysis_usage_examples:BadBinName', ...
            'Filename key-value pairs are invalid: %s', fname);
    end

    meta = struct();
    for i = 1:2:numel(kv)
        k = matlab.lang.makeValidName(kv{i});
        v = str2double(kv{i+1});
        if ~isnan(v)
            meta.(k) = v;
        end
    end

    keyOrder = {'dx','dy','dz','Lx','Ly','Lz'};
    vals = [];
    for i = 1:numel(keyOrder)
        k = keyOrder{i};
        if isfield(meta, k)
            vals(end+1) = meta.(k); %#ok<AGROW>
        end
    end

    if numel(vals) < 3
        error('analysis_usage_examples:NotEnoughVars', ...
            'Need at least 3 vars among dx/dy/dz/Lx/Ly/Lz in %s', fname);
    end

    dV = vals(1) * vals(2) * vals(3);
end
